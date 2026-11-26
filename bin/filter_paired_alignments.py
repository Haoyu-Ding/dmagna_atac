#!/usr/bin/env python3
"""Filter paired-end BAM alignments as complete fragments."""

from __future__ import annotations

import argparse
import csv
import itertools
import sys
from collections import Counter
from pathlib import Path

import pysam


def _edit_distance(read: pysam.AlignedSegment) -> int | None:
    try:
        return int(read.get_tag("NM"))
    except KeyError:
        return None


def _is_multimapped(read: pysam.AlignedSegment) -> bool:
    if read.mapping_quality < 1:
        return True
    try:
        return int(read.get_tag("NH")) > 1
    except KeyError:
        return False


def _is_soft_clipped(read: pysam.AlignedSegment) -> bool:
    return any(operation == 4 for operation, _length in (read.cigartuples or []))


def _is_fr_pair(read1: pysam.AlignedSegment, read2: pysam.AlignedSegment) -> bool:
    if read1.is_unmapped or read2.is_unmapped:
        return False
    if read1.reference_start < read2.reference_start:
        left, right = read1, read2
    elif read2.reference_start < read1.reference_start:
        left, right = read2, read1
    elif read1.template_length > 0:
        left, right = read1, read2
    elif read2.template_length > 0:
        left, right = read2, read1
    else:
        left, right = (read1, read2) if read1.is_read1 else (read2, read1)
    return not left.is_reverse and right.is_reverse


def _pair_failures(
    read1: pysam.AlignedSegment,
    read2: pysam.AlignedSegment,
    mitochondrial_contig: str,
    max_mismatches: int,
    max_insert_size: int,
) -> list[str]:
    failures: list[str] = []
    pair = (read1, read2)

    if any(read.is_unmapped or read.mate_is_unmapped for read in pair):
        failures.append("unmapped")
    if any(read.is_duplicate for read in pair):
        failures.append("duplicate")
    if any(read.is_qcfail for read in pair):
        failures.append("qc_fail")
    if any(_is_multimapped(read) for read in pair):
        failures.append("multi_mapped")
    if any(
        _edit_distance(read) is None or _edit_distance(read) > max_mismatches
        for read in pair
    ):
        failures.append("too_many_mismatches")
    if any(_is_soft_clipped(read) for read in pair):
        failures.append("soft_clipped")
    if max(abs(read1.template_length), abs(read2.template_length)) > max_insert_size:
        failures.append("insert_size_over_limit")
    if read1.reference_id != read2.reference_id:
        failures.append("different_chromosomes")
    if not _is_fr_pair(read1, read2):
        failures.append("not_FR_orientation")
    if any(
        read.reference_name == mitochondrial_contig
        or read.next_reference_name == mitochondrial_contig
        for read in pair
    ):
        failures.append("mitochondrial")
    return failures


def filter_bam(
    input_bam: Path,
    output_bam: Path,
    metrics_path: Path,
    mitochondrial_contig: str,
    max_mismatches: int,
    max_insert_size: int,
) -> None:
    counts: Counter[str] = Counter()
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    metrics_path.parent.mkdir(parents=True, exist_ok=True)

    with pysam.AlignmentFile(input_bam, "rb") as source, pysam.AlignmentFile(
        output_bam, "wb", template=source
    ) as destination:
        alignments = source.fetch(until_eof=True)
        for _query_name, group in itertools.groupby(alignments, key=lambda read: read.query_name):
            records = list(group)
            counts["pairs_seen"] += 1
            counts["non_primary_alignments_removed"] += sum(
                read.is_secondary or read.is_supplementary for read in records
            )

            primary = [
                read for read in records if not read.is_secondary and not read.is_supplementary
            ]
            read1s = [read for read in primary if read.is_read1]
            read2s = [read for read in primary if read.is_read2]
            if len(read1s) != 1 or len(read2s) != 1:
                counts["pairs_incomplete_or_ambiguous"] += 1
                counts["pairs_removed"] += 1
                continue

            first, second = read1s[0], read2s[0]
            failures = _pair_failures(
                first,
                second,
                mitochondrial_contig,
                max_mismatches,
                max_insert_size,
            )
            for reason in failures:
                counts[f"pairs_failed_{reason}"] += 1
            if failures:
                counts["pairs_removed"] += 1
                continue

            destination.write(first)
            destination.write(second)
            counts["pairs_retained"] += 1

    with metrics_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        for name in sorted(counts):
            writer.writerow([name, counts[name]])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--metrics", required=True, type=Path)
    parser.add_argument("--mitochondrial-contig", required=True)
    parser.add_argument("--max-mismatches", type=int, default=4)
    parser.add_argument("--max-insert-size", type=int, default=2000)
    args = parser.parse_args()

    try:
        filter_bam(
            args.input,
            args.output,
            args.metrics,
            args.mitochondrial_contig,
            args.max_mismatches,
            args.max_insert_size,
        )
    except (OSError, ValueError) as exc:
        print(f"ERROR: paired-alignment filtering failed: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
