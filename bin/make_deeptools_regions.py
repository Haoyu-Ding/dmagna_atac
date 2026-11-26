#!/usr/bin/env python3
"""Create gene-level BED regions from a matching GFF3 or GTF annotation."""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path
from typing import TextIO

import pysam


def _open_annotation(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _attributes(value: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    if "=" in value:
        for entry in value.strip().strip(";").split(";"):
            if "=" in entry:
                key, content = entry.split("=", 1)
                attributes[key] = content
    else:
        for key, content in re.findall(r'([^;\s]+)\s+"([^"]*)"', value):
            attributes[key] = content
    return attributes


def _feature_rows(annotation: Path):
    with _open_annotation(annotation) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            seqid, _source, feature, start, end, _score, strand, _phase, attrs = fields
            try:
                start0 = int(start) - 1
                end0 = int(end)
            except ValueError:
                continue
            if start0 < 0 or end0 <= start0 or strand not in {"+", "-"}:
                continue
            yield seqid, feature.lower(), start0, end0, strand, _attributes(attrs)


def _gene_regions(annotation: Path) -> list[tuple[str, int, int, str, str]]:
    gene_features: list[tuple[str, int, int, str, str]] = []
    transcripts: dict[str, tuple[str, int, int, str, str]] = {}

    for seqid, feature, start, end, strand, attrs in _feature_rows(annotation):
        if feature == "gene":
            identifier = attrs.get("Name") or attrs.get("gene_id") or attrs.get("ID")
            identifier = identifier or f"{seqid}:{start}-{end}"
            gene_features.append((seqid, start, end, strand, identifier))
        elif feature in {"mrna", "transcript", "ncrna", "rrna", "trna"}:
            gene_id = attrs.get("gene_id")
            parent = attrs.get("Parent")
            if not gene_id and parent:
                gene_id = parent.split(",")[0]
            identifier = gene_id or attrs.get("ID") or attrs.get("transcript_id")
            identifier = identifier or f"{seqid}:{start}-{end}"
            name = attrs.get("gene_name") or attrs.get("Name") or identifier
            prior = transcripts.get(identifier)
            if prior is None:
                transcripts[identifier] = (seqid, start, end, strand, name)
            else:
                if prior[0] != seqid or prior[3] != strand:
                    raise ValueError(f"conflicting transcript locations for {identifier}")
                transcripts[identifier] = (
                    seqid,
                    min(prior[1], start),
                    max(prior[2], end),
                    strand,
                    prior[4],
                )

    if gene_features:
        return gene_features
    return list(transcripts.values())


def write_regions(annotation: Path, bam: Path, output: Path) -> None:
    regions = _gene_regions(annotation)
    if not regions:
        raise ValueError("annotation contains no gene, mRNA, or transcript features")

    with pysam.AlignmentFile(bam, "rb") as alignments:
        reference_lengths = dict(zip(alignments.references, alignments.lengths))
        reference_names = set(reference_lengths)
    matched = [region for region in regions if region[0] in reference_names]
    unmatched = len(regions) - len(matched)
    if not matched:
        raise ValueError(
            "none of the annotation sequence names match the filtered BAM reference; "
            "use annotation built for the exact reference FASTA"
        )
    if unmatched:
        unmatched_names = sorted({region[0] for region in regions if region[0] not in reference_names})
        preview = ", ".join(unmatched_names[:5])
        raise ValueError(
            f"{unmatched} annotated gene regions use sequence names absent from the BAM "
            f"({preview}); use annotation built for the exact reference FASTA"
        )
    out_of_bounds = [
        (seqid, start, end, reference_lengths[seqid])
        for seqid, start, end, _strand, _name in matched
        if end > reference_lengths[seqid]
    ]
    if out_of_bounds:
        seqid, start, end, contig_length = out_of_bounds[0]
        raise ValueError(
            f"annotation coordinate {seqid}:{start}-{end} exceeds BAM reference length "
            f"{contig_length}; use annotation built for the exact reference FASTA"
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        for seqid, start, end, strand, name in matched:
            safe_name = re.sub(r"\s+", "_", name)
            handle.write(f"{seqid}\t{start}\t{end}\t{safe_name}\t0\t{strand}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation", required=True, type=Path)
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    try:
        write_regions(args.annotation, args.bam, args.output)
    except (OSError, ValueError) as exc:
        print(f"ERROR: cannot create deepTools gene regions: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
