#!/usr/bin/env python3
import argparse, csv, re, shutil, subprocess, sys
from pathlib import Path
from collections import defaultdict

try:
    from .logging_utils import print_log_only, run_logged_command
except ImportError:
    from logging_utils import print_log_only, run_logged_command

try:
    from .putative_blast import FilterSettings, run_putative_blastp
except ImportError:
    from putative_blast import FilterSettings, run_putative_blastp

# This script runs hmmscan on all the HMMs with input protein-coding files and outputs detailed TSV 
# btexgenie_summary with columns: sample, hmm, hits, scores, evalues, cutoff_used, total_genes, hit_headers.
# The cutoff_used column is is the sequence-level GA threshold from the HMM library.
# The evalue filter is applied to the sequence-level evalue.
#

DNA_EXTENSIONS = {".fna", ".fa", ".fasta"}
PROTEIN_EXTENSIONS = {".faa", ".aa", ".pep"}
DOMTBL_SUBDIR = "hmmscan_output"


def check_bin(name):
    from shutil import which
    if which(name) is None:
        raise SystemExit(f"[err] '{name}' not found in PATH")

def read_names_from_hmm_lib(hmm_lib: Path) -> list[str]:
    """
    Parse concatenated HMM library to extract model names (lines starting with 'NAME').
    """
    names = []
    if not hmm_lib.exists():
        return names
    with open(hmm_lib, "rt", errors="ignore") as fh:
        for line in fh:
            if line.startswith("NAME"):
                parts = line.strip().split()
                if len(parts) >= 2:
                    names.append(parts[1])
    return names


def read_ga_thresholds_from_hmm_lib(hmm_lib: Path) -> dict[str, tuple[float, float]]:
    """
    Parse model-specific GA thresholds from a concatenated HMM library.
    Returns: model_name -> (seq_ga_bits, dom_ga_bits)
    """
    ga_by_model: dict[str, tuple[float, float]] = {}
    if not hmm_lib.exists():
        return ga_by_model

    cur_name: str | None = None
    with open(hmm_lib, "rt", errors="ignore") as fh:
        for line in fh:
            if line.startswith("NAME"):
                parts = line.strip().split()
                cur_name = parts[1] if len(parts) >= 2 else None
                continue
            if cur_name and line.startswith("GA"):
                nums = re.findall(r"[-+]?\d*\.?\d+", line)
                if not nums:
                    continue
                if len(nums) == 1:
                    seq_ga = dom_ga = float(nums[0])
                else:
                    seq_ga = float(nums[0])
                    dom_ga = float(nums[1])
                ga_by_model[cur_name] = (seq_ga, dom_ga)
    return ga_by_model


def filter_domtbl_to_ga(
    in_domtbl: Path,
    out_domtbl: Path,
    ga_by_model: dict[str, tuple[float, float]],
) -> None:
    """
    Create GA-filtered domtblout from an unfiltered domtblout using model GA
    thresholds parsed from the HMM library.
    """
    warned_missing_ga: set[str] = set()
    with open(in_domtbl, "rt") as fr, open(out_domtbl, "wt") as fw:
        for line in fr:
            if line.startswith("#") or not line.strip():
                fw.write(line)
                continue
            parts = re.split(r"\s+", line.strip())
            if len(parts) < 22:
                continue
            hmm_name = parts[0]
            try:
                full_bits = float(parts[7])   # full sequence score
                dom_bits = float(parts[13])   # this domain score
            except ValueError:
                continue
            ga = ga_by_model.get(hmm_name)
            if ga is None:
                if hmm_name not in warned_missing_ga:
                    print(f"[warn] no GA threshold found in HMM library for model '{hmm_name}'; keeping row")
                    warned_missing_ga.add(hmm_name)
                fw.write(line)
                continue
            seq_ga, dom_ga = ga
            if full_bits >= seq_ga and dom_bits >= dom_ga:
                fw.write(line)


def parse_domtbl(path: Path):
    """
    Returns list of (seqid, hmm_name, seq_e, bits) tuples from an HMMER hmmscan domtblout.

    For hmmscan domain table:
      target name (col 1)   -> HMM (model) name
      query name  (col 4)   -> sequence id
      full seq E-value (col 7) -> parts[6]
      full seq score   (col 8) -> parts[7]
    """
    rows = []
    with open(path, "rt") as fh:
        for ln, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = re.split(r"\s+", line.strip())
            if len(parts) < 22:
                continue
            hmm_name = parts[0]     # target name = HMM
            seqid = parts[3]        # query name  = protein ID
            seq_e = float(parts[6]) # full sequence E-value
            bits = float(parts[7])  # full sequence bit score
            rows.append((seqid, hmm_name, seq_e, bits))
    return rows


def find_input_files(input_path: Path) -> list[Path]:
    if input_path.is_file():
        return [input_path]
    if input_path.is_dir():
        return [p for p in sorted(input_path.iterdir()) if p.is_file()]
    return []


def validate_protein_fasta(path: Path) -> tuple[bool, str]:
    """
    Validate that a file is FASTA and contains only protein-like residue symbols.
    Returns (ok, reason_if_not_ok).
    """
    allowed = set("ABCDEFGHIKLMNPQRSTUVWXYZ*-.")
    nuc_chars = set("ACGTUN")
    saw_header = False
    saw_seq = False
    total_alpha = 0
    nuc_alpha = 0
    saw_protein_specific = False
    try:
        with open(path, "rt", errors="ignore") as fh:
            for raw in fh:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    saw_header = True
                    continue
                if not saw_header:
                    return False, "missing FASTA header line starting with '>'"
                saw_seq = True
                seq = line.replace(" ", "").upper()
                bad = [c for c in seq if c not in allowed]
                if bad:
                    return False, f"contains non-protein residue character(s): {''.join(sorted(set(bad)))}"
                for c in seq:
                    if c.isalpha():
                        total_alpha += 1
                        if c in nuc_chars:
                            nuc_alpha += 1
                        else:
                            saw_protein_specific = True
    except OSError as e:
        return False, f"cannot read file: {e}"

    if not saw_header:
        return False, "no FASTA headers found"
    if not saw_seq:
        return False, "no sequence lines found"
    # Reject likely nucleotide FASTA: sequence alphabet is almost entirely A/C/G/T/U/N
    # and we never observe protein-specific residue letters.
    if total_alpha > 0:
        nuc_frac = nuc_alpha / total_alpha
        if nuc_frac >= 0.95 and not saw_protein_specific:
            return False, "sequence composition looks nucleotide (A/C/G/T/U/N only)"
    return True, ""


def validate_dna_fasta(path: Path) -> tuple[bool, str]:
    """
    Validate that a file is FASTA and contains only nucleotide-like residue symbols.
    Returns (ok, reason_if_not_ok).
    """
    allowed = set("ACGTRYSWKMBDHVNU-.")
    saw_header = False
    saw_seq = False
    saw_alpha = False
    try:
        with open(path, "rt", errors="ignore") as fh:
            for raw in fh:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    saw_header = True
                    continue
                if not saw_header:
                    return False, "missing FASTA header line starting with '>'"
                saw_seq = True
                seq = line.replace(" ", "").upper()
                bad = [c for c in seq if c not in allowed]
                if bad:
                    return False, f"contains non-DNA residue character(s): {''.join(sorted(set(bad)))}"
                if any(c.isalpha() for c in seq):
                    saw_alpha = True
    except OSError as e:
        return False, f"cannot read file: {e}"

    if not saw_header:
        return False, "no FASTA headers found"
    if not saw_seq:
        return False, "no sequence lines found"
    if not saw_alpha:
        return False, "no nucleotide sequence characters found"
    return True, ""


def detect_sequence_kind(path: Path) -> tuple[str | None, str]:
    suffix = path.suffix.lower()
    if suffix in PROTEIN_EXTENSIONS:
        ok, reason = validate_protein_fasta(path)
        return ("protein", "") if ok else (None, reason)
    if suffix == ".fna":
        ok, reason = validate_dna_fasta(path)
        return ("dna", "") if ok else (None, reason)
    if suffix in DNA_EXTENSIONS:
        protein_ok, _ = validate_protein_fasta(path)
        if protein_ok:
            return "protein", ""
        dna_ok, dna_reason = validate_dna_fasta(path)
        if dna_ok:
            return "dna", ""
        return None, dna_reason
    return None, f"unsupported file extension: {suffix or '<none>'}"


def classify_inputs(input_files: list[Path]) -> tuple[str, list[Path]]:
    classified: list[tuple[Path, str]] = []
    skipped = 0
    for path in input_files:
        kind, reason = detect_sequence_kind(path)
        if kind is None:
            print(f"[warn] skipping {path.name}: {reason}")
            skipped += 1
            continue
        classified.append((path, kind))

    if not classified:
        raise SystemExit("[err] no valid FASTA inputs found after validation")

    kinds = {kind for _, kind in classified}
    if len(kinds) != 1:
        raise SystemExit("[err] mixed DNA and protein inputs in one run are not supported")

    input_kind = classified[0][1]
    kept_files = [path for path, _ in classified]
    if skipped:
        print(f"[info] retained {len(kept_files)} validated {input_kind} input file(s)")
    return input_kind, kept_files


def run_prodigal_for_inputs(
    genome_fastas: list[Path],
    output_root: Path,
    prodigal_mode: str,
) -> tuple[list[Path], Path]:
    prodigal_root = output_root / "prodigal_output"
    prodigal_root.mkdir(parents=True, exist_ok=True)

    protein_fastas: list[Path] = []
    for genome_fna in genome_fastas:
        sample = genome_fna.stem
        sample_dir = prodigal_root / sample
        sample_dir.mkdir(parents=True, exist_ok=True)
        out_faa = sample_dir / f"{sample}.faa"
        out_gbk = sample_dir / f"{sample}_prodigal.gbk"

        print(f"[*] Prodigal gene calling ({prodigal_mode}) on {genome_fna.name}")
        cmd = [
            "prodigal",
            "-i",
            str(genome_fna),
            "-a",
            str(out_faa),
            "-o",
            str(out_gbk),
            "-p",
            prodigal_mode,
        ]
        run_logged_command(cmd)
        protein_fastas.append(out_faa)

    return protein_fastas, prodigal_root


def materialize_protein_inputs(
    protein_inputs: list[Path],
    output_root: Path,
) -> tuple[list[Path], Path]:
    prodigal_root = output_root / "prodigal_output"
    prodigal_root.mkdir(parents=True, exist_ok=True)

    protein_fastas: list[Path] = []
    for protein_src in protein_inputs:
        sample = protein_src.stem
        sample_dir = prodigal_root / sample
        sample_dir.mkdir(parents=True, exist_ok=True)
        out_faa = sample_dir / f"{sample}.faa"
        shutil.copyfile(protein_src, out_faa)
        protein_fastas.append(out_faa)

    return protein_fastas, prodigal_root


def load_headers_for_fasta(protein_faa: Path) -> dict[str, str]:
    """
    Build mapping from sequence id -> full header line (including leading '>') for a single FASTA.
    The id is the first whitespace-delimited token after '>'.
    """
    mp: dict[str, str] = {}
    try:
        with open(protein_faa, "rt", errors="ignore") as fh:
            for ln in fh:
                if not ln.startswith(">"):
                    continue
                hdr = ln.rstrip("\n\r")
                seqid = hdr[1:].split()[0]
                mp.setdefault(seqid, hdr)  # keep first occurrence
    except FileNotFoundError:
        pass
    return mp


def run_hmmscan_for_sample(
    protein_faa: Path,
    hmm_lib: Path,
    domtbl_dir: Path,
    cpus: int,
    ga_by_model: dict[str, tuple[float, float]],
):
    """
    Run hmmscan once per sample against the concatenated HMM library.
    Writes per-sample domtbl outputs:
      - <sample>.pre_ga.domtblout  (no --cut_ga)
      - <sample>.ga.domtblout      (filtered from pre-GA using model GA thresholds)
    """
    if not protein_faa or not protein_faa.exists():
        print(f"[warn] protein FASTA missing for {protein_faa} — skipping hmmscan")
        return

    if not hmm_lib.exists():
        raise SystemExit(f"[err] HMM library not found: {hmm_lib}")

    domtbl_dir.mkdir(parents=True, exist_ok=True)
    sample = protein_faa.stem
    out_dom_pre = domtbl_dir / f"{sample}.pre_ga.domtblout"
    out_dom_ga = domtbl_dir / f"{sample}.ga.domtblout"

    # 1) Always write unfiltered domtbl (before GA).
    print_log_only(f"[info] running hmmscan on {protein_faa.name}")
    cmd_pre = ["hmmscan", "--cpu", str(cpus)]
    cmd_pre += ["--domtblout", str(out_dom_pre), str(hmm_lib), str(protein_faa)]
    try:
        run_logged_command(cmd_pre)
    except subprocess.CalledProcessError as e:
        print(
            f"[err] hmmscan (pre-GA) failed for {protein_faa.name} with code {e.returncode}",
            file=sys.stderr,
        )
        return

    # 2) Derive GA-filtered domtbl from pre-GA output (no second hmmscan run).
    filter_domtbl_to_ga(out_dom_pre, out_dom_ga, ga_by_model)


def run_kofam_for_inputs(proteins_path: Path, cpus: int) -> None:
    kofam_script = Path(__file__).resolve().with_name("kofamscan.py")
    if not kofam_script.exists():
        raise SystemExit(f"[err] KOfam wrapper not found: {kofam_script}")
    if not proteins_path.is_dir():
        print(
            f"[warn] skipping KOfam search because input is not a directory: {proteins_path}"
        )
        return

    cmd = [
        sys.executable,
        str(kofam_script),
        "--genomes-dir",
        str(proteins_path.resolve()),
        "--cpus",
        str(cpus),
    ]
    print(f"[*] starting KOfam search on {proteins_path}")
    run_logged_command(cmd)


def main(argv: list[str] | None = None):
    ap = argparse.ArgumentParser(
        description=(
            "Run hmmscan vs protein FASTAs using a concatenated HMM library, "
            "writing one row per best protein hit as "
            "sample,hmm,score,evalue,cutoff_used,hit_header. "
            "Sample name is taken from the FASTA basename. "
            "GA thresholds are read from each HMM in the concatenated library."
        )
    )
    ap.add_argument("--hmm-lib", required=True, help="concatenated HMM library for hmmscan")
    ap.add_argument(
        "--genomes-dir",
        dest="genomes_dir",
        required=True,
        help=(
            "Directory or single file containing either genome DNA FASTA files or protein FASTA files. "
            "Input files must all be the same type."
        ),
    )
    ap.add_argument("--cpus", type=int, default=8, help="threads for hmmscan")
    ap.add_argument(
        "--evalue",
        type=float,
        default=1e-5,
        help="full-sequence E-value cutoff applied to output hits after GA filtering (default: 1e-5)",
    )
    prodigal_group = ap.add_mutually_exclusive_group()
    prodigal_group.add_argument(
        "-meta",
        dest="prodigal_mode",
        action="store_const",
        const="meta",
        help="use Prodigal meta mode when DNA genomes are provided",
    )
    prodigal_group.add_argument(
        "-single",
        dest="prodigal_mode",
        action="store_const",
        const="single",
        help="use Prodigal single mode when DNA genomes are provided",
    )
    ap.add_argument(
        "--kofam",
        action="store_true",
        help="Run the KOfam search after hmmscan completes",
    )
    ap.add_argument(
        "--out",
        required=True,
        help="CSV path for hit-level output (sample,hmm,score,evalue,cutoff_used,hit_header)",
    )

    ap.add_argument(
        "--putative-blast",
        action="store_true",
        help="Search proposed BTEX pathway proteins against the input protein sequences",
    )

    ap.add_argument(
        "--putative-targets",
        help="FASTA file or directory containing putative target proteins",
    )

    ap.add_argument(
        "--putative-metadata",
        help="Metadata TSV describing the putative target proteins",
    )

    ap.add_argument(
        "--blast-min-identity",
        type=float,
        default=25.0,
    )

    ap.add_argument(
        "--blast-min-query-coverage",
        type=float,
        default=50.0,
        help=(
            "Minimum reference-protein coverage percent for filtered putative-target "
            "BLASTP hits. Reference proteins are the BLASTP queries (default: 50.0)"
        ),
    )

    ap.add_argument(
        "--blast-evalue",
        type=float,
        default=1e-5,
        help=(
            "BLASTP E-value threshold for putative target searches "
            "(default: 1e-5)"
        ),
    )
    args = ap.parse_args(argv)

    check_bin("hmmscan")

    hmm_lib = Path(args.hmm_lib)
    genomes_dir = Path(args.genomes_dir)
    ga_by_model = read_ga_thresholds_from_hmm_lib(hmm_lib)

    out_path = Path(args.out)
    output_root = out_path.parent if out_path.parent != Path("") else Path(".")

    if not hmm_lib.exists():
        raise SystemExit(f"[err] HMM library not found: {hmm_lib}")
    all_model_names = read_names_from_hmm_lib(hmm_lib)
    if not all_model_names:
        raise SystemExit(f"[err] no model names found in HMM library: {hmm_lib}")
    model_set = set(all_model_names)

    input_files = find_input_files(genomes_dir)
    if not input_files:
        raise SystemExit(
            f"[err] no input files found in {genomes_dir}"
        )
    input_kind, validated_inputs = classify_inputs(input_files)
    if input_kind == "dna":
        if args.prodigal_mode is None:
            raise SystemExit(
                "[err] DNA genome inputs detected. Use either --meta or --single to choose the Prodigal mode."
            )
        check_bin("prodigal")
        protein_fastas, kofam_input_dir = run_prodigal_for_inputs(
            genome_fastas=validated_inputs,
            output_root=output_root,
            prodigal_mode=args.prodigal_mode,
        )
    else:
        protein_fastas, kofam_input_dir = materialize_protein_inputs(
            protein_inputs=validated_inputs,
            output_root=output_root,
        )
        if args.prodigal_mode is not None:
            print_log_only("[warn] ignoring Prodigal mode flag because protein FASTA inputs were provided")

    header_maps: dict[str, dict[str, str]] = {}
    samples: list[str] = []

    # Phase 1: ensure per-sample domtbl exists (one hmmscan per sample)
    for protein_faa in protein_fastas:
        sample = protein_faa.stem
        samples.append(sample)

        domtbl_dir = output_root / DOMTBL_SUBDIR / sample
        run_hmmscan_for_sample(
            protein_faa=protein_faa,
            hmm_lib=hmm_lib,
            domtbl_dir=domtbl_dir,
            cpus=args.cpus,
            ga_by_model=ga_by_model,
        )

        mp = load_headers_for_fasta(protein_faa)
        header_maps[sample] = mp

    # Phase 2: parse domtbls and build one output row per best protein hit
    output_rows: list[list[str | float]] = []
    count_rows: dict[tuple[str, str], dict[str, str | float | int]] = {}
    warned_missing_ga = set()

    for sample in samples:
        domtbl_dir = output_root / DOMTBL_SUBDIR / sample
        domtbl_path = domtbl_dir / f"{sample}.ga.domtblout"
        if not domtbl_path.exists():
            print(f"[info] no domtbl for {sample}; skipping")
            continue

        seqid_to_header = header_maps.get(sample, {})

        hits_per_protein = defaultdict(list)
        for seqid, model_name, seq_e, bits in parse_domtbl(domtbl_path):
            if seq_e > args.evalue:
                continue
            model = model_name
            if model not in model_set:
                continue

            ga = ga_by_model.get(model)
            ga_seq = ga[0] if ga is not None else None
            if ga_seq is None and model not in warned_missing_ga:
                print(
                    f"[warn] no GA threshold found in HMM library for model '{model}'; "
                    "using only other filters for this model"
                )
                warned_missing_ga.add(model)

            diff_score = bits - (ga_seq if ga_seq is not None else 0.0)
            hits_per_protein[seqid].append(
                {
                    "model": model,
                    "seqid": seqid,
                    "diff": diff_score,
                    "bits": bits,
                    "seq_e": seq_e,
                    "ga_seq": ga_seq,
                }
            )

        for seqid, hit_list in hits_per_protein.items():
            if not hit_list:
                continue
            best_hit = max(hit_list, key=lambda h: h["diff"])
            hdr = seqid_to_header.get(seqid, f">{seqid}")
            output_rows.append(
                [
                    sample,
                    best_hit["model"],
                    best_hit["bits"],
                    best_hit["seq_e"],
                    best_hit["ga_seq"] if best_hit["ga_seq"] is not None else "",
                    hdr,
                ]
            )
            key = (sample, best_hit["model"])
            if key not in count_rows:
                count_rows[key] = {
                    "sample": sample,
                    "hmm": best_hit["model"],
                    "score": str(best_hit["bits"]),
                    "evalue": str(best_hit["seq_e"]),
                    "cutoff_used": best_hit["ga_seq"] if best_hit["ga_seq"] is not None else "",
                    "count": 0,
                }
            count_rows[key]["count"] = int(count_rows[key]["count"]) + 1

    # Phase 3: write output
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as fw:
        w = csv.writer(fw)
        w.writerow(
            [
                "sample",
                "hmm",
                "score",
                "evalue",
                "cutoff_used",
                "hit_header",
            ]
        )
        for row in sorted(output_rows, key=lambda r: (str(r[0]), str(r[1]), str(r[5]))):
            w.writerow(row)

    print(f"wrote {out_path}")

    counts_path = out_path.with_name("btex_genie_summary_counts.csv")
    with open(counts_path, "w", newline="") as fw:
        w = csv.writer(fw)
        w.writerow(
            [
                "sample",
                "hmm",
                "score",
                "evalue",
                "cutoff_used",
                "count",
            ]
        )
        for key in sorted(count_rows):
            row = count_rows[key]
            w.writerow(
                [
                    row["sample"],
                    row["hmm"],
                    row["score"],
                    row["evalue"],
                    row["cutoff_used"],
                    row["count"],
                ]
            )

    print(f"wrote {counts_path}")

    if args.putative_blast:
        if not args.putative_targets:
            raise SystemExit(
                "[err] --putative-targets is required when "
                "--putative-blast is enabled"
            )

        putative_filters = FilterSettings(
            evalue=args.blast_evalue,
            min_identity=args.blast_min_identity,
            min_query_coverage=args.blast_min_query_coverage,
        )

        try:
            run_putative_blastp(
                protein_fastas=protein_fastas,
                targets_path=Path(args.putative_targets),
                output_root=output_root,
                cpus=args.cpus,
                metadata_path=(
                    Path(args.putative_metadata)
                    if args.putative_metadata
                    else None
                ),
                filters=putative_filters,
            )
        except (FileNotFoundError, ValueError, RuntimeError) as exc:
            raise SystemExit(f"[err] putative BLASTP search failed: {exc}") from exc
    else:
        print(
            "[info] skipping putative-target BLASTP; "
            "use --putative-blast to enable it"
        )

    if args.kofam:
        try:
            run_kofam_for_inputs(kofam_input_dir, args.cpus)
        except subprocess.CalledProcessError as e:
            raise SystemExit(
                f"[err] KOfam search failed with exit code {e.returncode}"
            ) from e
    else:
        print("[info] skipping KOfam search by default; use --kofam to enable it")


if __name__ == "__main__":
    main()
