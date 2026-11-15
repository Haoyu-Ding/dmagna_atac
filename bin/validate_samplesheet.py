#!/usr/bin/env python3
"""Validate and normalize the explicit input manifest for the ATAC-seq pipeline."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


REQUIRED_COLUMNS = ("sample_id", "fastq_1", "fastq_2")
OPTIONAL_COLUMNS = ("condition", "replicate")
FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz")
SAMPLE_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
READ_TOKEN_PATTERN = re.compile(r"(?:^|[_\-.])R([12])(?=$|[_\-.])", re.IGNORECASE)


class ValidationResult:
    def __init__(self, rows: List[dict], warnings: List[str], errors: List[str]):
        self.rows = rows
        self.warnings = warnings
        self.errors = errors


def _clean(value: Optional[str]) -> str:
    return (value or "").strip()


def _resolve_path(value: str, samplesheet_dir: Path) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = samplesheet_dir / path
    return path.resolve()


def _fastq_signature(path: Path) -> Tuple[int, str]:
    """Return read number and filename with the R1/R2 token removed."""
    lower_name = path.name.lower()
    suffix = next((item for item in FASTQ_SUFFIXES if lower_name.endswith(item)), None)
    if suffix is None:
        raise ValueError("filename does not end in .fastq.gz or .fq.gz")

    stem = path.name[: -len(suffix)]
    matches = list(READ_TOKEN_PATTERN.finditer(stem))
    if len(matches) != 1:
        raise ValueError("filename must contain exactly one R1/R2 token")

    match = matches[0]
    read_number = int(match.group(1))
    sequencing_unit = (stem[: match.start()] + stem[match.end() :]).lower()
    return read_number, sequencing_unit


def _check_gzip(path: Path) -> bool:
    try:
        with path.open("rb") as handle:
            return handle.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def validate_samplesheet(input_path: Path) -> ValidationResult:
    warnings: List[str] = []
    errors: List[str] = []
    rows: List[dict] = []
    samplesheet_dir = input_path.parent

    try:
        handle = input_path.open("r", encoding="utf-8-sig", newline="")
    except OSError as exc:
        return ValidationResult([], [], [f"cannot open sample sheet {input_path}: {exc}"])

    with handle:
        reader = csv.DictReader(handle)
        header = reader.fieldnames
        if not header:
            return ValidationResult([], [], ["sample sheet is empty or has no header row"])
        if any(not _clean(column) for column in header):
            errors.append("sample sheet contains an empty column name")
        if len(header) != len(set(header)):
            errors.append("sample sheet contains duplicate column names")

        missing = [column for column in REQUIRED_COLUMNS if column not in header]
        if missing:
            errors.append("missing required column(s): " + ", ".join(missing))
            return ValidationResult([], warnings, errors)

        unknown = [
            column
            for column in header
            if column not in REQUIRED_COLUMNS and column not in OPTIONAL_COLUMNS
        ]
        if unknown:
            warnings.append(
                "unknown column(s) will be recorded but ignored by workflow logic: "
                + ", ".join(unknown)
            )

        metadata_by_sample: Dict[str, Tuple[str, str, int]] = {}
        sample_pairs = defaultdict(set)
        sample_r1_to_r2 = defaultdict(dict)
        sample_r2_to_r1 = defaultdict(dict)
        used_files: Dict[Path, Tuple[int, str, str]] = {}

        for line_number, raw_row in enumerate(reader, start=2):
            if raw_row is None or not any(_clean(value) for value in raw_row.values()):
                continue
            if None in raw_row:
                errors.append(f"line {line_number}: contains more fields than the header")
                continue

            sample_id = _clean(raw_row.get("sample_id"))
            fastq_1_value = _clean(raw_row.get("fastq_1"))
            fastq_2_value = _clean(raw_row.get("fastq_2"))
            condition = _clean(raw_row.get("condition"))
            replicate = _clean(raw_row.get("replicate"))

            row_errors: List[str] = []
            if not sample_id:
                row_errors.append("sample_id is empty")
            elif not SAMPLE_ID_PATTERN.fullmatch(sample_id):
                row_errors.append(
                    "sample_id must match ^[A-Za-z0-9][A-Za-z0-9._-]*$ and be safe for output filenames"
                )
            if not fastq_1_value or not fastq_2_value:
                row_errors.append("both fastq_1 and fastq_2 are required")

            fastq_1 = _resolve_path(fastq_1_value, samplesheet_dir) if fastq_1_value else None
            fastq_2 = _resolve_path(fastq_2_value, samplesheet_dir) if fastq_2_value else None
            for label, path in (("fastq_1", fastq_1), ("fastq_2", fastq_2)):
                if path is None:
                    continue
                if path.suffixes[-2:] not in [[".fastq", ".gz"], [".fq", ".gz"]]:
                    row_errors.append(f"{label} must end in .fastq.gz or .fq.gz: {path}")
                elif not path.is_file():
                    row_errors.append(f"{label} does not exist or is not a file: {path}")
                elif not _check_gzip(path):
                    row_errors.append(f"{label} is not gzip-compressed: {path}")

            if fastq_1 is not None and fastq_2 is not None:
                try:
                    read_1, unit_1 = _fastq_signature(fastq_1)
                    read_2, unit_2 = _fastq_signature(fastq_2)
                    if read_1 != 1 or read_2 != 2:
                        row_errors.append("fastq_1/fastq_2 must be an R1/R2 pair in that order")
                    elif unit_1 != unit_2:
                        row_errors.append(
                            "fastq_1 and fastq_2 do not belong to the same sequencing unit"
                        )
                except ValueError as exc:
                    row_errors.append(f"cannot validate R1/R2 pairing: {exc}")

            if row_errors:
                errors.extend(f"line {line_number}: {message}" for message in row_errors)
                continue

            assert fastq_1 is not None and fastq_2 is not None
            pair = (fastq_1, fastq_2)
            if pair in sample_pairs[sample_id]:
                errors.append(f"line {line_number}: duplicated read pair for sample_id {sample_id}")
            sample_pairs[sample_id].add(pair)

            previous_r2 = sample_r1_to_r2[sample_id].get(fastq_1)
            if previous_r2 is not None and previous_r2 != fastq_2:
                errors.append(f"line {line_number}: sample_id {sample_id} reuses R1 with a different R2")
            sample_r1_to_r2[sample_id][fastq_1] = fastq_2

            previous_r1 = sample_r2_to_r1[sample_id].get(fastq_2)
            if previous_r1 is not None and previous_r1 != fastq_1:
                errors.append(f"line {line_number}: sample_id {sample_id} reuses R2 with a different R1")
            sample_r2_to_r1[sample_id][fastq_2] = fastq_1

            for path, role in ((fastq_1, "fastq_1"), (fastq_2, "fastq_2")):
                previous = used_files.get(path)
                if previous is not None:
                    previous_line, previous_sample, previous_role = previous
                    errors.append(
                        f"line {line_number}: {role} {path} is already used on line "
                        f"{previous_line} as {previous_role} for sample_id {previous_sample}"
                    )
                else:
                    used_files[path] = (line_number, sample_id, role)

            metadata = (condition, replicate)
            previous_metadata = metadata_by_sample.get(sample_id)
            if previous_metadata is not None and previous_metadata[:2] != metadata:
                errors.append(
                    f"line {line_number}: repeated sample_id {sample_id} has inconsistent "
                    f"condition/replicate metadata (first seen on line {previous_metadata[2]})"
                )
            else:
                metadata_by_sample.setdefault(sample_id, (condition, replicate, line_number))

            normalized_row = {column: _clean(raw_row.get(column)) for column in header}
            normalized_row["sample_id"] = sample_id
            normalized_row["fastq_1"] = str(fastq_1)
            normalized_row["fastq_2"] = str(fastq_2)
            normalized_row["condition"] = condition if "condition" in header else ""
            normalized_row["replicate"] = replicate if "replicate" in header else ""
            rows.append(normalized_row)

    if not rows:
        errors.append("no valid sample rows remain after validation")

    return ValidationResult(rows, warnings, errors)


def _write_outputs(
    input_path: Path,
    result: ValidationResult,
    validated_output: Path,
    input_files_output: Path,
    report_output: Path,
) -> None:
    validated_output.parent.mkdir(parents=True, exist_ok=True)
    input_files_output.parent.mkdir(parents=True, exist_ok=True)
    report_output.parent.mkdir(parents=True, exist_ok=True)

    with input_path.open("r", encoding="utf-8-sig", newline="") as handle:
        header = csv.DictReader(handle).fieldnames or list(REQUIRED_COLUMNS)

    if result.rows:
        with validated_output.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=header, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(result.rows)

    lane_numbers = defaultdict(int)
    with input_files_output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "lane", "fastq_1", "fastq_2", "condition", "replicate"])
        for row in result.rows:
            lane_numbers[row["sample_id"]] += 1
            writer.writerow(
                [
                    row["sample_id"],
                    lane_numbers[row["sample_id"]],
                    row["fastq_1"],
                    row["fastq_2"],
                    row.get("condition", ""),
                    row.get("replicate", ""),
                ]
            )

    status = "PASS" if not result.errors else "FAIL"
    with report_output.open("w", encoding="utf-8") as handle:
        handle.write(f"status: {status}\n")
        handle.write(f"input: {input_path.resolve()}\n")
        handle.write(f"valid_rows: {len(result.rows)}\n")
        handle.write(f"unique_samples: {len({row['sample_id'] for row in result.rows})}\n")
        handle.write("\n")
        if result.warnings:
            handle.write("warnings:\n")
            handle.writelines(f"- {warning}\n" for warning in result.warnings)
            handle.write("\n")
        if result.errors:
            handle.write("errors:\n")
            handle.writelines(f"- {error}\n" for error in result.errors)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--validated-output", type=Path, required=True)
    parser.add_argument("--input-files-output", type=Path, required=True)
    parser.add_argument("--report-output", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    result = validate_samplesheet(args.input)
    _write_outputs(
        args.input,
        result,
        args.validated_output,
        args.input_files_output,
        args.report_output,
    )
    if result.errors:
        print("Sample sheet validation failed:", file=sys.stderr)
        for error in result.errors:
            print(f"- {error}", file=sys.stderr)
        return 2
    print(f"Validated {len(result.rows)} read-pair row(s) across "
          f"{len({row['sample_id'] for row in result.rows})} sample(s).")
    for warning in result.warnings:
        print(f"Warning: {warning}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
