import csv
import gzip
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "bin" / "validate_samplesheet.py"
SPEC = importlib.util.spec_from_file_location("validate_samplesheet", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class ValidateSamplesheetTests(unittest.TestCase):
    def make_fastq_pair(self, directory: Path, stem: str):
        r1 = directory / f"{stem}_R1_001.fastq.gz"
        r2 = directory / f"{stem}_R2_001.fastq.gz"
        for path in (r1, r2):
            with gzip.open(path, "wt", encoding="ascii") as handle:
                handle.write("@read1\nACGT\n+\n!!!!\n")
        return r1, r2

    def write_sheet(self, directory: Path, rows, extra_header=False):
        sheet = directory / "samplesheet.csv"
        header = ["sample_id", "fastq_1", "fastq_2", "condition", "replicate"]
        if extra_header:
            header.append("notes")
        with sheet.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=header)
            writer.writeheader()
            writer.writerows(rows)
        return sheet

    def test_valid_rows_are_normalized_and_groupable(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            r1, r2 = self.make_fastq_pair(directory, "sampleA_L001")
            r1b, r2b = self.make_fastq_pair(directory, "sampleA_L002")
            sheet = self.write_sheet(
                directory,
                [
                    {"sample_id": "sample_A", "fastq_1": r1.name, "fastq_2": r2.name, "condition": "fish", "replicate": "1"},
                    {"sample_id": "sample_A", "fastq_1": r1b.name, "fastq_2": r2b.name, "condition": "fish", "replicate": "1"},
                ],
            )
            result = MODULE.validate_samplesheet(sheet)
            self.assertEqual(result.errors, [])
            self.assertEqual(len(result.rows), 2)
            self.assertEqual(result.rows[0]["fastq_1"], str(r1.resolve()))

    def test_mismatched_pair_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            r1, _ = self.make_fastq_pair(directory, "sampleA_L001")
            _, r2 = self.make_fastq_pair(directory, "sampleB_L001")
            sheet = self.write_sheet(
                directory,
                [{"sample_id": "sample_A", "fastq_1": r1.name, "fastq_2": r2.name, "condition": "", "replicate": ""}],
            )
            result = MODULE.validate_samplesheet(sheet)
            self.assertTrue(any("same sequencing unit" in error for error in result.errors))

    def test_inconsistent_metadata_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            r1, r2 = self.make_fastq_pair(directory, "sampleA_L001")
            r1b, r2b = self.make_fastq_pair(directory, "sampleA_L002")
            sheet = self.write_sheet(
                directory,
                [
                    {"sample_id": "sample_A", "fastq_1": r1.name, "fastq_2": r2.name, "condition": "fish", "replicate": "1"},
                    {"sample_id": "sample_A", "fastq_1": r1b.name, "fastq_2": r2b.name, "condition": "control", "replicate": "1"},
                ],
            )
            result = MODULE.validate_samplesheet(sheet)
            self.assertTrue(any("inconsistent" in error for error in result.errors))

    def test_unknown_columns_are_warned_about(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            r1, r2 = self.make_fastq_pair(directory, "sampleA_L001")
            sheet = self.write_sheet(
                directory,
                [{"sample_id": "sample_A", "fastq_1": r1.name, "fastq_2": r2.name, "condition": "fish", "replicate": "1", "notes": "keep"}],
                extra_header=True,
            )
            result = MODULE.validate_samplesheet(sheet)
            self.assertEqual(result.errors, [])
            self.assertTrue(any("unknown column" in warning for warning in result.warnings))


if __name__ == "__main__":
    unittest.main()
