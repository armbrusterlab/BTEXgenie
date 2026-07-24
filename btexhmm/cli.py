import argparse
import shlex
from pathlib import Path

from .hmmscan import main as hmmscan_main
from .logging_utils import command_logger

HERE = Path(__file__).resolve().parent
DEFAULT_HMM_LIB = HERE / "hmms" / "BTEXgenie_updated.hmm"
DEFAULT_PUTATIVE_TARGETS = HERE / "data" / "putative_targets.faa"
DEFAULT_PUTATIVE_METADATA = HERE / "data" / "putative_targets.tsv"


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Annotate a directory or single genome DNA FASTA files or protein FASTA files "
            "with the BTEX HMM database. DNA inputs are first gene-called with Prodigal."
        )
    )
    p.add_argument(
        "-g",
        "--genome-dir",
        dest="genome_dir",
        metavar="GENOMES",
        required=True,
        help=(
            "Directory or single file containing genome DNA FASTA or protein FASTA input "
            "(for example *.fna, *.fa, *.fasta, *.faa)."
        ),
    )
    p.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="Output directory for results",
    )
    p.add_argument(
        "--cpus",
        type=int,
        default=8,
        help="Number of CPUs for hmmscan (default: 8)",
    )
    prodigal_group = p.add_mutually_exclusive_group()
    prodigal_group.add_argument(
        "--meta",
        dest="prodigal_mode",
        action="store_const",
        const="meta",
        help="Use Prodigal meta mode for DNA genome inputs",
    )
    prodigal_group.add_argument(
        "--single",
        dest="prodigal_mode",
        action="store_const",
        const="single",
        help="Use Prodigal single mode for DNA genome inputs (default mode)",
    )
    p.add_argument(
        "--evalue",
        type=float,
        default=1e-5,
        help="Full-sequence E-value cutoff applied to output hits after GA filtering (default: 1e-5)",
    )
    p.add_argument(
        "--kofam",
        action="store_true",
        help="Run the HMM search against the KOfam database in addition to the BTEX-HMM database.",
    )

    p.add_argument(
        "--blast-min-identity",
        type=float,
        default=25.0,
        help="Minimum percent identity for reporting putative-target BLASTP hits (default: 25.0)",
    )
    p.add_argument(
        "--blast-min-query-coverage",
        type=float,
        default=50.0,
        help=(
            "Minimum reference-protein coverage percent for filtered putative-target "
            "BLASTP hits. Reference proteins are the BLASTP queries (default: 50.0)"
        ),
    )
    p.add_argument(
        "--blast-evalue",
        type=float,
        default=1e-5,
        help=(
            "BLASTP E-value threshold for putative target searches "
            "(default: 1e-5)"
        ),
    )
    return p.parse_args()


def main():
    args = parse_args()
    prodigal_mode = args.prodigal_mode or "single"

    if not DEFAULT_HMM_LIB.is_file():
        raise SystemExit(
            f"[err] bundled BTEXgenie HMM database not found: {DEFAULT_HMM_LIB}"
        )

    genomes = Path(args.genome_dir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    out_csv = outdir / "btex_genie_summary.csv"
    hmmscan_argv = [
        "--hmm-lib", str(DEFAULT_HMM_LIB),
        "--genomes-dir", str(genomes),
        "--out", str(out_csv),
        "--cpus", str(args.cpus),
        "--putative-blast",
        "--putative-targets",
        str(DEFAULT_PUTATIVE_TARGETS),
        "--putative-metadata",
        str(DEFAULT_PUTATIVE_METADATA),
        "--blast-min-identity",
        str(args.blast_min_identity),
        "--blast-min-query-coverage",
        str(args.blast_min_query_coverage),
        "--blast-evalue",
        str(args.blast_evalue),
    ]
    if prodigal_mode == "meta":
        hmmscan_argv.append("-meta")
    elif prodigal_mode == "single":
        hmmscan_argv.append("-single")
    if args.evalue is not None:
        hmmscan_argv.extend(["--evalue", str(args.evalue)])
    if args.kofam:
        hmmscan_argv.append("--kofam")

    top_cmd = [
        "btex-annotate",
        "-g",
        str(genomes),
        "-o",
        str(outdir),
        "--cpus",
        str(args.cpus),
    ]
    if prodigal_mode == "meta":
        top_cmd.append("--meta")
    elif prodigal_mode == "single":
        top_cmd.append("--single")
    if args.evalue is not None:
        top_cmd.extend(["--evalue", str(args.evalue)])
    if args.kofam:
        top_cmd.append("--kofam")
    if args.blast_min_identity != 25.0:
        top_cmd.extend(["--blast-min-identity", str(args.blast_min_identity)])
    if args.blast_min_query_coverage != 50.0:
        top_cmd.extend(["--blast-min-query-coverage", str(args.blast_min_query_coverage)])
    if args.blast_evalue != 1e-5:
        top_cmd.extend(["--blast-evalue", str(args.blast_evalue)])

    log_path = outdir / "btex_annotate.log"
    with command_logger(log_path):
        print(f"[info] writing btex-annotate log to {log_path}")
        print(f"[cmd] {shlex.join(top_cmd)}")
        print(f"[info] using BTEXgenie HMM database: {DEFAULT_HMM_LIB}")
        print("[info] running BTEXgenie hmmscan and putative-target BLASTP on input genomes")
        hmmscan_main(hmmscan_argv)
        print("Done!")
