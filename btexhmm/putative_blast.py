#!/usr/bin/env python3
"""BLASTP search module for putative BTEX pathway proteins.

This module is designed to be imported by the BTEXgenie annotation workflow after
protein FASTA files have been generated or materialized. It can also be run as a
standalone command for testing.

The putative reference proteins are used as BLASTP queries and each user proteome
is used as the subject. The module writes only two final output files:

1. btex_putative_blastp_filtered_summary.csv
   One row per passing hit with compact reference-aware columns.

2. btex_putative_blastp_unfiltered_summary.csv
   One row per BLASTP hit before applying the user identity and coverage
   thresholds. These rows are still limited by the BLASTP search E-value used
   in the command.

"""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence

try:
    from .logging_utils import print_log_only, run_logged_command
except ImportError:
    try:
        from logging_utils import print_log_only, run_logged_command
    except ImportError:
        def print_log_only(message: str, end: str = "\n") -> None:
            print(message, end=end)

        def run_logged_command(cmd: Sequence[str | Path]) -> None:
            printable = " ".join(str(part) for part in cmd)
            print(f"[cmd] {printable}")
            subprocess.run([str(part) for part in cmd], check=True)


PROTEIN_EXTENSIONS = {".faa", ".fa", ".fasta", ".aa", ".pep"}

BLAST_FIELDS = [
    "qseqid",
    "sseqid",
    "pident",
    "ppos",
    "nident",
    "positive",
    "length",
    "mismatch",
    "gapopen",
    "qlen",
    "slen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]

SUMMARY_FIELDS = [
    "sample",
    "evalue",
    "percent_identity",
    "percent_similarity",
    "reference_coverage",
    "reference_subject_length",
    "bitscore",
    "query_hit_header",
    "reference_hit_header",
]


@dataclass(frozen=True)
class TargetMetadata:
    query_id: str
    system: str
    component: str
    report_name: str
    evidence_level: str


@dataclass(frozen=True)
class FilterSettings:
    evalue: float = 1e-5
    min_identity: float = 25.0
    min_query_coverage: float = 50.0
    min_subject_coverage: float = 50.0


@dataclass(frozen=True)
class FastaRecord:
    record_id: str
    header: str
    sequence: str


def check_bin(name: str) -> None:
    if shutil.which(name) is None:
        raise SystemExit(f"[err] '{name}' not found in PATH")


def iter_fasta_records(path: Path) -> Iterator[FastaRecord]:
    """Yield FASTA records using the first header token as the record ID."""
    current_header: str | None = None
    sequence_parts: list[str] = []

    with open(path, "rt", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    record_id = current_header.split()[0]
                    sequence = "".join(sequence_parts).replace(" ", "").upper()
                    if not sequence:
                        raise ValueError(f"Empty FASTA sequence in {path}: {record_id}")
                    yield FastaRecord(record_id, current_header, sequence)
                current_header = line[1:].strip()
                if not current_header:
                    raise ValueError(f"Blank FASTA header in {path}")
                sequence_parts = []
            else:
                if current_header is None:
                    raise ValueError(
                        f"Sequence encountered before a FASTA header in {path}"
                    )
                sequence_parts.append(line)

    if current_header is not None:
        record_id = current_header.split()[0]
        sequence = "".join(sequence_parts).replace(" ", "").upper()
        if not sequence:
            raise ValueError(f"Empty FASTA sequence in {path}: {record_id}")
        yield FastaRecord(record_id, current_header, sequence)


def find_protein_fastas(path: Path) -> list[Path]:
    """Return protein FASTA files from a single file or a flat directory."""
    if path.is_file():
        return [path]
    if not path.is_dir():
        return []
    return [
        item
        for item in sorted(path.iterdir())
        if item.is_file() and item.suffix.lower() in PROTEIN_EXTENSIONS
    ]


def build_query_fasta(targets_path: Path, output_fasta: Path) -> list[str]:
    """Create one query FASTA from a target file or directory.

    For a target directory, each FASTA file must contain exactly one sequence and
    the file stem is used as the query ID. This keeps output identifiers stable
    even when source FASTA headers contain long database descriptions.

    For a single multi-FASTA input, the original first-token record IDs are used.
    """
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    query_ids: list[str] = []
    seen_ids: set[str] = set()

    if targets_path.is_file():
        source_files = [targets_path]
        use_file_stem = False
    elif targets_path.is_dir():
        source_files = [
            item
            for item in sorted(targets_path.iterdir())
            if item.is_file() and item.suffix.lower() in PROTEIN_EXTENSIONS
        ]
        use_file_stem = True
    else:
        raise FileNotFoundError(f"Putative target path not found: {targets_path}")

    if not source_files:
        raise ValueError(f"No protein FASTA files found in {targets_path}")

    with open(output_fasta, "wt", encoding="utf-8") as out_handle:
        for source in source_files:
            records = list(iter_fasta_records(source))
            if use_file_stem and len(records) != 1:
                raise ValueError(
                    f"Expected exactly one FASTA record in {source}, found {len(records)}"
                )

            for record in records:
                query_id = source.stem if use_file_stem else record.record_id
                if query_id in seen_ids:
                    raise ValueError(f"Duplicate putative target ID: {query_id}")
                seen_ids.add(query_id)
                query_ids.append(query_id)

                original_description = record.header
                out_handle.write(f">{query_id} source_header={original_description}\n")
                for start in range(0, len(record.sequence), 80):
                    out_handle.write(record.sequence[start : start + 80] + "\n")

    return query_ids


def load_target_metadata(
    metadata_path: Path | None,
    query_ids: Iterable[str],
) -> dict[str, TargetMetadata]:
    """Load optional target metadata and provide safe defaults for missing rows."""
    query_id_set = set(query_ids)
    metadata: dict[str, TargetMetadata] = {}

    if metadata_path is not None:
        if not metadata_path.exists():
            raise FileNotFoundError(f"Putative target metadata not found: {metadata_path}")

        with open(metadata_path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"query_id", "system", "component", "report_name", "evidence_level"}
            missing = required.difference(reader.fieldnames or [])
            if missing:
                raise ValueError(
                    "Metadata TSV is missing required column(s): "
                    + ", ".join(sorted(missing))
                )

            for row in reader:
                query_id = (row.get("query_id") or "").strip()
                if not query_id:
                    continue
                if query_id in metadata:
                    raise ValueError(f"Duplicate metadata row for query ID: {query_id}")
                metadata[query_id] = TargetMetadata(
                    query_id=query_id,
                    system=(row.get("system") or query_id).strip(),
                    component=(row.get("component") or query_id).strip(),
                    report_name=(row.get("report_name") or query_id).strip(),
                    evidence_level=(row.get("evidence_level") or "putative").strip(),
                )

        unknown_metadata_ids = set(metadata).difference(query_id_set)
        if unknown_metadata_ids:
            print_log_only(
                "[warn] metadata contains query IDs not present in the target FASTA: "
                + ", ".join(sorted(unknown_metadata_ids))
            )

    for query_id in query_id_set:
        metadata.setdefault(
            query_id,
            TargetMetadata(
                query_id=query_id,
                system=query_id,
                component=query_id,
                report_name=f"{query_id}-like protein",
                evidence_level="putative",
            ),
        )

    return metadata


def load_fasta_headers(path: Path) -> dict[str, str]:
    headers: dict[str, str] = {}
    with open(path, "rt", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            full_header = line.rstrip("\r\n")
            record_id = full_header[1:].split()[0]
            headers.setdefault(record_id, full_header)
    return headers


def run_blastp_for_sample(
    protein_faa: Path,
    query_faa: Path,
    temp_dir: Path,
    cpus: int,
    evalue: float,
) -> Path:
    """Run BLASTP for one sample and return the temporary TSV path."""
    if not protein_faa.exists():
        raise FileNotFoundError(f"Protein FASTA not found: {protein_faa}")
    if not query_faa.exists():
        raise FileNotFoundError(f"Putative query FASTA not found: {query_faa}")

    sample = protein_faa.stem
    raw_path = temp_dir / f"{sample}.putative_blastp.tsv"

    cmd = [
        "blastp",
        "-query",
        str(query_faa),
        "-subject",
        str(protein_faa),
        "-evalue",
        str(evalue),
        "-seg",
        "yes",
        "-soft_masking",
        "true",
        "-max_hsps",
        "1",
        "-num_threads",
        str(max(1, cpus)),
        "-outfmt",
        "6 " + " ".join(BLAST_FIELDS),
        "-out",
        str(raw_path),
    ]

    print_log_only(f"[info] running putative-target BLASTP on {protein_faa.name}")
    try:
        run_logged_command(cmd)
    except subprocess.CalledProcessError as exc:
        raw_path.unlink(missing_ok=True)
        raise RuntimeError(
            f"BLASTP failed for {protein_faa.name} with exit code {exc.returncode}"
        ) from exc

    return raw_path


def _to_float(value: str, field_name: str, source: Path) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid {field_name} value in {source}: {value!r}") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"Non-finite {field_name} value in {source}: {value!r}")
    return parsed


def _to_int(value: str, field_name: str, source: Path) -> int:
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid {field_name} value in {source}: {value!r}") from exc


def passes_user_thresholds(
    row: dict[str, str | float | int],
    filters: FilterSettings,
) -> bool:
    return (
        float(row["evalue"]) <= filters.evalue
        and float(row["percent_identity"]) >= filters.min_identity
        and float(row["reference_coverage"]) >= filters.min_query_coverage
    )


def parse_blastp_hits(
    raw_path: Path,
    sample: str,
    metadata: dict[str, TargetMetadata],
    sample_headers: dict[str, str],
    reference_headers: dict[str, str],
) -> list[dict[str, str | float | int]]:
    """Parse a headerless BLASTP TSV and calculate reference-protein coverage."""
    rows: list[dict[str, str | float | int]] = []

    with open(raw_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, fieldnames=BLAST_FIELDS, delimiter="\t")
        for blast_row in reader:
            if not blast_row or not blast_row.get("qseqid"):
                continue
            missing_values = [field for field in BLAST_FIELDS if blast_row.get(field) is None]
            if missing_values:
                raise ValueError(
                    f"Raw BLASTP file {raw_path} has a malformed row missing: "
                    + ", ".join(missing_values)
                )

            query_id = str(blast_row["qseqid"])
            subject_id = str(blast_row["sseqid"])
            target_meta = metadata.get(
                query_id,
                TargetMetadata(
                    query_id=query_id,
                    system=query_id,
                    component=query_id,
                    report_name=f"{query_id}-like protein",
                    evidence_level="putative",
                ),
            )

            pident = _to_float(str(blast_row["pident"]), "pident", raw_path)
            ppos = _to_float(str(blast_row["ppos"]), "ppos", raw_path)
            reference_length = _to_int(str(blast_row["qlen"]), "qlen", raw_path)
            qstart = _to_int(str(blast_row["qstart"]), "qstart", raw_path)
            qend = _to_int(str(blast_row["qend"]), "qend", raw_path)
            evalue = _to_float(str(blast_row["evalue"]), "evalue", raw_path)
            bitscore = _to_float(str(blast_row["bitscore"]), "bitscore", raw_path)

            reference_span = abs(qend - qstart) + 1
            reference_coverage = 100.0 * reference_span / reference_length

            rows.append(
                {
                    "sample": sample,
                    "system": target_meta.system,
                    "component": target_meta.component,
                    "query_id": query_id,
                    "subject_id": subject_id,
                    "report_name": target_meta.report_name,
                    "evidence_level": target_meta.evidence_level,
                    "evalue": evalue,
                    "percent_identity": round(pident, 3),
                    "percent_similarity": round(ppos, 3),
                    "reference_coverage": round(reference_coverage, 3),
                    "reference_subject_length": reference_length,
                    "bitscore": bitscore,
                    "query_hit_header": sample_headers.get(subject_id, f">{subject_id}"),
                    "reference_hit_header": reference_headers.get(query_id, f">{query_id}"),
                }
            )

    return rows


def write_csv_rows(
    path: Path,
    rows: Iterable[dict[str, object]],
    fieldnames: Sequence[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sort_hit_rows(rows: Iterable[dict[str, object]]) -> list[dict[str, object]]:
    return sorted(
        rows,
        key=lambda row: (
            str(row["sample"]),
            str(row["query_id"]),
            float(row["evalue"]),
            -float(row["bitscore"]),
            str(row["subject_id"]),
        ),
    )


def run_putative_blastp(
    protein_fastas: Sequence[Path],
    targets_path: Path,
    output_root: Path,
    cpus: int = 8,
    metadata_path: Path | None = None,
    filters: FilterSettings | None = None,
) -> dict[str, Path]:
    """Run the complete putative-target BLASTP workflow.

    This is the function intended to be imported by BTEXgenie's main annotation
    workflow after its ``protein_fastas`` list has been created.
    """
    if not protein_fastas:
        raise ValueError("No protein FASTA files were supplied to putative BLASTP")

    check_bin("blastp")
    filters = filters or FilterSettings()
    output_root.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict[str, str | float | int]] = []

    with tempfile.TemporaryDirectory(
        prefix=".btex_putative_blastp_tmp_",
        dir=output_root,
    ) as temp_name:
        temp_dir = Path(temp_name)
        query_faa = temp_dir / "putative_targets.combined.faa"
        query_ids = build_query_fasta(targets_path, query_faa)
        metadata = load_target_metadata(metadata_path, query_ids)
        reference_headers = load_fasta_headers(query_faa)

        for protein_faa in protein_fastas:
            sample = protein_faa.stem
            raw_path = run_blastp_for_sample(
                protein_faa=protein_faa,
                query_faa=query_faa,
                temp_dir=temp_dir,
                cpus=cpus,
                evalue=filters.evalue,
            )
            sample_headers = load_fasta_headers(protein_faa)
            all_rows.extend(
                parse_blastp_hits(
                    raw_path=raw_path,
                    sample=sample,
                    metadata=metadata,
                    sample_headers=sample_headers,
                    reference_headers=reference_headers,
                )
            )

    unfiltered_rows = sort_hit_rows(all_rows)
    filtered_rows = sort_hit_rows(
        row for row in all_rows if passes_user_thresholds(row, filters)
    )

    filtered_summary_path = output_root / "btex_putative_blastp_filtered_summary.csv"
    unfiltered_summary_path = output_root / "btex_putative_blastp_unfiltered_summary.csv"

    write_csv_rows(filtered_summary_path, filtered_rows, SUMMARY_FIELDS)
    write_csv_rows(unfiltered_summary_path, unfiltered_rows, SUMMARY_FIELDS)

    print(f"wrote {filtered_summary_path}")
    print(f"wrote {unfiltered_summary_path}")

    return {
        "filtered_summary": filtered_summary_path,
        "unfiltered_summary": unfiltered_summary_path,
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Search putative BTEX pathway proteins against one or more protein "
            "FASTAs with BLASTP and write filtered and unfiltered summary reports."
        )
    )
    parser.add_argument(
        "--proteins",
        required=True,
        help="Protein FASTA file or flat directory containing protein FASTAs",
    )
    parser.add_argument(
        "--targets",
        required=True,
        help=(
            "Putative target FASTA or directory of one-sequence-per-file target FASTAs"
        ),
    )
    parser.add_argument(
        "--metadata",
        help=(
            "Optional target metadata TSV with query_id, system, component, "
            "report_name, and evidence_level columns"
        ),
    )
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--cpus", type=int, default=8, help="BLASTP threads")
    parser.add_argument(
        "--evalue",
        type=float,
        default=1e-5,
        help=(
            "BLASTP search E-value threshold. The unfiltered summary includes all "
            "BLASTP rows returned under this search E-value before identity and "
            "coverage filtering."
        ),
    )
    parser.add_argument("--min-identity", type=float, default=25.0)
    parser.add_argument(
        "--min-query-coverage",
        type=float,
        default=50.0,
        help=(
            "Minimum reference protein coverage for filtered summary rows. "
            "The reference proteins are BLASTP queries in this implementation."
        ),
    )
    parser.add_argument(
        "--min-subject-coverage",
        type=float,
        default=50.0,
        help=(
            "Deprecated compatibility option. Sample subject coverage is no longer "
            "written to the summary reports or used for filtering."
        ),
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    protein_fastas = find_protein_fastas(Path(args.proteins))
    if not protein_fastas:
        raise SystemExit(f"[err] no protein FASTA files found in {args.proteins}")

    metadata_path = Path(args.metadata) if args.metadata else None
    filters = FilterSettings(
        evalue=args.evalue,
        min_identity=args.min_identity,
        min_query_coverage=args.min_query_coverage,
        min_subject_coverage=args.min_subject_coverage,
    )

    try:
        run_putative_blastp(
            protein_fastas=protein_fastas,
            targets_path=Path(args.targets),
            output_root=Path(args.outdir),
            cpus=args.cpus,
            metadata_path=metadata_path,
            filters=filters,
        )
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        raise SystemExit(f"[err] {exc}") from exc


if __name__ == "__main__":
    main()