"""Standard-library tests; all CSV fixtures and outputs are synthetic and temporary.

From the group repository root:
    py -B -m unittest discover -s tests -p "test_repair_solvent_pairs.py" -v
"""

from __future__ import annotations

import copy
import csv
import importlib.util
import io
import json
import shutil
import sys
import unittest
import uuid
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from pathlib import Path
from unittest.mock import patch


REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "python scripts" / "03_repair_solvent_pairs.py"
SPEC = importlib.util.spec_from_file_location("repair_solvent_pairs", SCRIPT)
repair = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = repair
SPEC.loader.exec_module(repair)

# Deliberately unordered and noncontiguous: IDs must not come from CSV row indices.
COMPOUNDS = [
    {"compound_id": "140", "smiles": "OC(C)C"},
    {"compound_id": "901", "smiles": "CC(=O)O"},
    {"compound_id": "29", "smiles": "O"},
    {"compound_id": "17", "smiles": "OC"},
]


def fixture(percent1=0, percent2=100, smiles1="O", smiles2="OC(C)C",
            target="-0.942900", temperature="283.15", label="arbitrary-label"):
    mapping = {row["smiles"]: int(row["compound_id"]) for row in COMPOUNDS}
    raw = {
        "Label": label, "LogS": target, "Temperature, K": temperature,
        "Drug": "Synthetic drug", "SMILES_0": "CC(=O)O", "SMILES_1": smiles1,
        "SMILES_2": smiles2, "comp_0": "100", "comp_1": str(percent1),
        "comp_2": str(percent2),
    }
    processed = {
        "value": target, "Temperature, K": temperature, "property": "Log solubility",
        "unit": "g/L", "cmp_ids_solvent": json.dumps(
            [mapping[s] for s in (smiles1, smiles2) if s in mapping]
        ),
        "cmp_mole_fractions_solvent": json.dumps(
            [float(x) / 100 for x in (percent1, percent2) if float(x) != 0]
        ),
        "cmp_ids_drug": "[901]", "cmp_mole_fractions_drug": "[1.0]",
        "note": 'Preserve this comma, quote " and newline\nexactly',
    }
    return raw, processed


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


@contextmanager
def dataset_fixture(records=None):
    # Keep all temporary file activity inside the authorized group repository.
    # Use inherited permissions: Python 3.13+ Windows mkdtemp's owner-only ACL
    # prevents a restricted sandbox token from accessing its own test fixtures.
    root = REPO_ROOT / f".repair-test-{uuid.uuid4().hex}"
    root.mkdir()
    try:
        records = records if records is not None else [fixture(), fixture(100, 0, "OC", "O")]
        paths = repair.input_paths(root)
        write_csv(paths["raw"], [raw for raw, processed in records])
        write_csv(paths["processed"], [processed for raw, processed in records])
        write_csv(paths["compounds"], COMPOUNDS)
        paths["metadata"].write_text('{"description": "Synthetic test dataset"}', encoding="utf-8")
        paths["upstream_processor"].write_text("# Synthetic sentinel; never executed\n", encoding="utf-8")
        yield root, paths
    finally:
        if root.resolve().parent != REPO_ROOT.resolve() or not root.name.startswith(".repair-test-"):
            raise ValueError("Refusing to clean up a fixture outside the group repository")
        shutil.rmtree(root)


def call_main(arguments):
    out, err = io.StringIO(), io.StringIO()
    with redirect_stdout(out), redirect_stderr(err):
        code = repair.main(arguments)
    return code, out.getvalue(), err.getvalue()


class RepairRowsTests(unittest.TestCase):
    def apply(self, raw, processed, compounds=None):
        return repair.repair_rows([raw], [processed], COMPOUNDS if compounds is None else compounds)

    def test_zero_solvent_1_preserves_solvent_2(self):
        raw, processed = fixture()
        result = self.apply(raw, processed)
        self.assertEqual(result.rows[0]["cmp_ids_solvent"], "[140]")
        self.assertEqual(result.changes[0]["zero_solvent_position"], 1)
        self.assertEqual(result.changes[0]["raw_label"], "arbitrary-label")

    def test_zero_solvent_2_preserves_solvent_1(self):
        raw, processed = fixture(100, 0, "OC", "O")
        result = self.apply(raw, processed)
        self.assertEqual(result.rows[0]["cmp_ids_solvent"], "[17]")
        self.assertEqual(result.changes[0]["zero_solvent_position"], 2)

    def test_only_id_cell_changes_and_inputs_are_not_mutated(self):
        raw, processed = fixture()
        processed["cmp_ids_solvent"] = "[29.0, 140.0]"
        processed["cmp_mole_fractions_solvent"] = "[1.000000]"
        processed["unit"] = " custom unit label "
        before = copy.deepcopy((raw, processed, COMPOUNDS))
        result = self.apply(raw, processed)
        self.assertEqual((raw, processed, COMPOUNDS), before)
        expected = processed.copy()
        expected["cmp_ids_solvent"] = "[140]"
        self.assertEqual(result.rows, [expected])

    def test_correct_mixture_is_preserved_exactly(self):
        raw, processed = fixture(25, 75)
        processed["cmp_ids_solvent"] = "[29.0, 140.0]"
        processed["cmp_mole_fractions_solvent"] = "[0.25000, 0.75000]"
        result = self.apply(raw, processed)
        self.assertEqual(result.rows, [processed])
        self.assertEqual(result.changes, [])

    def test_already_correct_explicit_zero_pairs_are_preserved(self):
        raw, processed = fixture()
        processed["cmp_mole_fractions_solvent"] = "[0.0, 1.0]"
        result = self.apply(raw, processed)
        self.assertEqual(result.rows, [processed])
        self.assertEqual(result.changes, [])

    def test_positive_small_fraction_is_never_dropped(self):
        raw, processed = fixture(1e-9, 100 - 1e-9)
        result = self.apply(raw, processed)
        self.assertEqual(result.rows, [processed])
        self.assertEqual(result.changes, [])

    def test_duplicates_record_count_and_order_are_preserved(self):
        a, b = fixture(), fixture(100, 0, "OC", "O", label="second")
        rows = [a, b, copy.deepcopy(a)]
        result = repair.repair_rows([r for r, p in rows], [p for r, p in rows], COMPOUNDS)
        self.assertEqual(len(result.rows), 3)
        self.assertEqual([r["cmp_ids_solvent"] for r in result.rows], ["[140]", "[17]", "[140]"])
        self.assertEqual(result.rows[0], result.rows[2])
        self.assertEqual([c["data_row"] for c in result.changes], [1, 2, 3])

    def test_second_pass_makes_no_changes(self):
        raw, processed = fixture()
        first = self.apply(raw, processed)
        second = repair.repair_rows([raw], first.rows, COMPOUNDS)
        self.assertEqual(second.rows, first.rows)
        self.assertEqual(second.changes, [])

    def test_duplicate_ids_or_smiles_are_ambiguous(self):
        raw, processed = fixture()
        additions = [
            {"compound_id": "29", "smiles": "new-smiles"},
            {"compound_id": "999", "smiles": "O"},
            COMPOUNDS[0].copy(),
        ]
        for addition in additions:
            with self.subTest(addition=addition), self.assertRaisesRegex(ValueError, "duplicate ID|ambiguous"):
                self.apply(raw, processed, COMPOUNDS + [addition])

    def test_invalid_mapping_ids_or_empty_smiles_are_rejected(self):
        raw, processed = fixture()
        for bad in ("-1", "1.5", "nan", "inf", "", "not-a-number"):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                self.apply(raw, processed, COMPOUNDS + [{"compound_id": bad, "smiles": "new"}])
        with self.assertRaisesRegex(ValueError, "empty SMILES"):
            self.apply(raw, processed, COMPOUNDS + [{"compound_id": "999", "smiles": " "}])

    def test_unmapped_positive_component_is_rejected(self):
        raw, processed = fixture()
        raw["SMILES_2"] = "unmapped"
        with self.assertRaisesRegex(ValueError, "positive raw component"):
            self.apply(raw, processed)

    def test_absent_zero_component_smiles_is_allowed_for_correct_row(self):
        raw, processed = fixture(smiles1="")
        result = self.apply(raw, processed)
        self.assertEqual(result.rows, [processed])
        self.assertEqual(result.changes, [])

    def test_unmapped_zero_component_cannot_justify_dropping_an_extra_id(self):
        raw, processed = fixture()
        for missing in ("", "unmapped"):
            with self.subTest(missing=missing):
                raw["SMILES_1"] = missing
                with self.assertRaisesRegex(ValueError, "solvent ID correspondence"):
                    self.apply(raw, processed)

    def test_unknown_fractional_or_repeated_processed_ids_are_rejected(self):
        raw, processed = fixture()
        for text in ("[29, 999]", "[29.5, 140]", "[29, 29]", "[True, 140]"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                processed["cmp_ids_solvent"] = text
                self.apply(raw, processed)

    def test_row_count_mismatch_is_rejected(self):
        raw, processed = fixture()
        with self.assertRaisesRegex(ValueError, "row count mismatch"):
            repair.repair_rows([raw], [processed, processed], COMPOUNDS)

    def test_swapped_rows_are_rejected(self):
        first = fixture(target="1.25")
        second = fixture(target="2.5")
        with self.assertRaisesRegex(ValueError, "correspondence failure for value"):
            repair.repair_rows([first[0], second[0]], [second[1], first[1]], COMPOUNDS)

    def test_temperature_misalignment_is_rejected(self):
        raw, processed = fixture()
        processed["Temperature, K"] = "298.15"
        with self.assertRaisesRegex(ValueError, "correspondence failure for Temperature"):
            self.apply(raw, processed)

    def test_drug_misalignment_is_rejected(self):
        raw, processed = fixture()
        processed["cmp_ids_drug"] = "[17]"
        with self.assertRaisesRegex(ValueError, "drug correspondence"):
            self.apply(raw, processed)

    def test_solvent_misalignment_with_identical_target_and_temperature_is_rejected(self):
        raw, processed = fixture()
        processed["cmp_ids_solvent"] = "[17, 140]"
        with self.assertRaisesRegex(ValueError, "solvent ID correspondence"):
            self.apply(raw, processed)

    def test_reversed_ids_are_rejected_before_repair(self):
        raw, processed = fixture()
        processed["cmp_ids_solvent"] = "[140, 29]"
        with self.assertRaisesRegex(ValueError, "solvent ID correspondence"):
            self.apply(raw, processed)

    def test_unexpected_fraction_difference_is_not_rounded_or_normalized(self):
        raw, processed = fixture()
        # Within the sum tolerance, but outside the much tighter pair tolerance.
        processed["cmp_mole_fractions_solvent"] = "[0.999999999]"
        with self.assertRaisesRegex(ValueError, "existing fraction list differs"):
            self.apply(raw, processed)

    def test_same_length_wrong_pairs_are_rejected(self):
        raw, processed = fixture(25, 75)
        processed["cmp_mole_fractions_solvent"] = "[0.75, 0.25]"
        with self.assertRaisesRegex(ValueError, "solvent pairs do not match"):
            self.apply(raw, processed)

    def test_mismatch_without_exactly_one_zero_is_rejected(self):
        raw, processed = fixture(25, 75)
        processed["cmp_mole_fractions_solvent"] = "[1.0]"
        with self.assertRaisesRegex(ValueError, "does not fit the zero-fraction defect"):
            self.apply(raw, processed)

    def test_other_length_mismatches_are_rejected(self):
        raw, processed = fixture()
        processed["cmp_ids_solvent"] = "[29, 140, 17]"
        with self.assertRaisesRegex(ValueError, "does not fit"):
            self.apply(raw, processed)

    def test_invalid_raw_compositions_are_rejected(self):
        raw, processed = fixture()
        for first, second in ((-1, 101), (50, 49), (0, 0), ("nan", 100)):
            with self.subTest(first=first, second=second), self.assertRaises(ValueError):
                raw["comp_1"], raw["comp_2"] = str(first), str(second)
                self.apply(raw, processed)
        raw, processed = fixture()
        raw["comp_0"] = "50"
        with self.assertRaisesRegex(ValueError, "invalid raw composition"):
            self.apply(raw, processed)

    def test_malformed_nonfinite_or_boolean_fraction_lists_are_rejected(self):
        raw, processed = fixture()
        for text in ("not a list", "(1.0,)", "[]", "[True]", "[1e309]", "[-1]", "[0.5]"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                processed["cmp_mole_fractions_solvent"] = text
                self.apply(raw, processed)

    def test_invalid_target_temperature_property_and_unit_are_rejected(self):
        raw, processed = fixture()
        for key, value in (("value", "nan"), ("Temperature, K", "0"),
                           ("property", "Other property"), ("unit", " ")):
            with self.subTest(key=key), self.assertRaises(ValueError):
                row = processed.copy()
                row[key] = value
                self.apply(raw, row)


class FileAndCliTests(unittest.TestCase):
    def test_default_and_explicit_dry_runs_write_nothing_and_label_validation(self):
        with dataset_fixture() as (root, paths):
            initial_files = set(root.rglob("*"))
            before = repair.fingerprints(list(paths.values()))
            for extra in ([], ["--dry-run"]):
                with self.subTest(extra=extra):
                    code, output, error = call_main(["--dataset-dir", str(root), *extra])
                    self.assertEqual(code, 0, error)
                    self.assertIn("ORIGINAL DATA VALIDATION (unchanged disk input)", output)
                    self.assertIn("SOLVENT_LENGTH_MISMATCH: 2", output)
                    self.assertIn("CORRECTED DATA VALIDATION (in-memory records", output)
                    self.assertIn("Structural errors: 0", output)
                    self.assertIn("Affected records detected from contents: 2", output)
                    self.assertIn("Second pass changes: 0", output)
                    self.assertIn("No repaired dataset or change log saved.", output)
                    self.assertEqual(set(root.rglob("*")), initial_files)
                    self.assertEqual(repair.fingerprints(list(paths.values())), before)

    def test_existing_validator_namespace_is_not_patched(self):
        with dataset_fixture() as (root, paths):
            namespace = repair.load_validator(paths["validator"])
            loader = namespace["validate"].__globals__["load_csv"]
            raw = repair.load_table(paths["raw"], repair.RAW_REQUIRED)
            processed = repair.load_table(paths["processed"], repair.PROCESSED_REQUIRED)
            fixed = repair.repair_rows(raw.rows, processed.rows, COMPOUNDS)
            _, corrected, _ = repair.validate_with_existing(namespace, paths, fixed.rows)
            _, original, _ = repair.validate_with_existing(namespace, paths)
            self.assertFalse(any(i["severity"] == "ERROR" for i in corrected))
            self.assertEqual(sum(i["severity"] == "ERROR" for i in original), 2)
            self.assertIs(namespace["validate"].__globals__["load_csv"], loader)

    def test_explicit_save_round_trip_and_second_saved_run(self):
        with dataset_fixture([fixture(), fixture()]) as (root, paths):
            before = repair.fingerprints(list(paths.values()))
            original = repair.load_table(paths["processed"], repair.PROCESSED_REQUIRED)
            output_csv, log = root / "outputs" / "repaired.csv", root / "outputs" / "changes.csv"
            code, summary, error = call_main([
                "--dataset-dir", str(root), "--output-csv", str(output_csv), "--change-log", str(log)
            ])
            self.assertEqual(code, 0, error)
            saved = repair.load_table(output_csv, repair.PROCESSED_REQUIRED)
            self.assertEqual(saved.header, original.header)
            self.assertEqual(len(saved.rows), 2)
            self.assertEqual(saved.rows[0], saved.rows[1])
            for source, corrected in zip(original.rows, saved.rows):
                expected = source.copy()
                expected["cmp_ids_solvent"] = "[140]"
                self.assertEqual(corrected, expected)
            with log.open(encoding="utf-8", newline="") as handle:
                changes = list(csv.DictReader(handle))
            self.assertEqual(len(changes), 2)
            self.assertEqual([c["data_row"] for c in changes], ["1", "2"])
            self.assertEqual(changes[0]["changed_column"], "cmp_ids_solvent")
            self.assertEqual(changes[0]["before"], original.rows[0]["cmp_ids_solvent"])
            self.assertEqual(changes[0]["after"], "[140]")
            self.assertEqual(json.loads(changes[0]["input_sha256"])[str(paths["raw"].resolve())],
                             before[str(paths["raw"].resolve())])
            self.assertIn("Repeated-condition candidate rows retained: 2", summary)
            output2, log2 = root / "second.csv", root / "second-changes.csv"
            code, summary, error = call_main([
                "--dataset-dir", str(root), "--processed-input", str(output_csv),
                "--output-csv", str(output2), "--change-log", str(log2)
            ])
            self.assertEqual(code, 0, error)
            self.assertIn("Affected records detected from contents: 0", summary)
            self.assertEqual(repair.load_table(output2, repair.PROCESSED_REQUIRED).rows, saved.rows)
            with log2.open(encoding="utf-8", newline="") as handle:
                self.assertEqual(list(csv.DictReader(handle)), [])
            self.assertEqual(repair.fingerprints(list(paths.values())), before)

    def test_output_arguments_must_be_paired_and_cannot_accompany_dry_run(self):
        for arguments in (["--output-csv", "unused.csv"], ["--change-log", "unused.csv"],
                          ["--dry-run", "--output-csv", "unused.csv", "--change-log", "unused-log.csv"]):
            with self.subTest(arguments=arguments), redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as stopped:
                    repair.main(arguments)
                self.assertEqual(stopped.exception.code, 2)

    def test_inputs_and_scripts_are_protected_in_both_output_positions(self):
        with dataset_fixture() as (root, paths):
            before = repair.fingerprints(list(paths.values()) + [SCRIPT])
            for protected in [*paths.values(), SCRIPT]:
                for csv_path, log_path in ((protected, root / "new-log.csv"),
                                           (root / "new.csv", protected)):
                    with self.subTest(protected=protected, csv_path=csv_path):
                        code, output, error = call_main([
                            "--dataset-dir", str(root), "--output-csv", str(csv_path),
                            "--change-log", str(log_path)
                        ])
                        self.assertEqual(code, 2)
                        self.assertIn("protected input", error)
                        self.assertFalse((root / "new.csv").exists())
                        self.assertFalse((root / "new-log.csv").exists())
            self.assertEqual(repair.fingerprints(list(paths.values()) + [SCRIPT]), before)

    def test_existing_outputs_are_protected_and_no_other_output_is_created(self):
        with dataset_fixture() as (root, paths):
            existing = root / "existing.csv"
            existing.write_bytes(b"existing output sentinel")
            for output_csv, log in ((existing, root / "new-dir" / "log.csv"),
                                     (root / "new-dir" / "fixed.csv", existing)):
                with self.subTest(output_csv=output_csv):
                    code, _, error = call_main([
                        "--dataset-dir", str(root), "--output-csv", str(output_csv), "--change-log", str(log)
                    ])
                    self.assertEqual(code, 2)
                    self.assertIn("already exists", error)
                    self.assertEqual(existing.read_bytes(), b"existing output sentinel")
                    self.assertFalse((root / "new-dir").exists())

    def test_same_destination_and_resolved_input_alias_are_rejected(self):
        with dataset_fixture() as (root, paths):
            with self.assertRaisesRegex(ValueError, "different paths"):
                repair.guard_outputs(root / "out.csv", root / "out.csv", list(paths.values()))
            alias = root / "raw_data" / ".." / "raw_data" / "DrugSolubilityData.csv"
            with self.assertRaisesRegex(ValueError, "protected input"):
                repair.guard_outputs(alias, root / "log.csv", list(paths.values()))

    def test_unexpected_error_in_later_row_prevents_any_saved_output(self):
        first, second = fixture(), fixture(100, 0, "OC", "O")
        second[1]["cmp_mole_fractions_solvent"] = "[0.999999999]"
        with dataset_fixture([first, second]) as (root, paths):
            before = repair.fingerprints(list(paths.values()))
            code, _, error = call_main([
                "--dataset-dir", str(root), "--output-csv", str(root / "new" / "fixed.csv"),
                "--change-log", str(root / "new" / "log.csv")
            ])
            self.assertEqual(code, 2)
            self.assertIn("data row 2", error)
            self.assertIn("existing fraction list differs", error)
            self.assertFalse((root / "new").exists())
            self.assertEqual(repair.fingerprints(list(paths.values())), before)

    def test_unrelated_original_validation_error_stops_before_correction(self):
        raw, processed = fixture()
        processed["value"] = "123"
        with dataset_fixture([(raw, processed)]) as (root, paths):
            code, output, error = call_main(["--dataset-dir", str(root), "--dry-run"])
            self.assertEqual(code, 2)
            self.assertIn("TARGET_NOT_PRESERVED", error)
            self.assertNotIn("CORRECTED DATA VALIDATION", output)

    def test_write_failure_removes_only_new_partial_outputs(self):
        raw, processed = fixture()
        result = repair.repair_rows([raw], [processed], COMPOUNDS)
        with dataset_fixture() as (root, paths):
            before = repair.fingerprints(list(paths.values()))
            csv_path, log = root / "fixed.csv", root / "changes.csv"
            with patch.object(repair.csv.DictWriter, "writerows", side_effect=OSError("disk full")):
                with self.assertRaisesRegex(OSError, "disk full"):
                    repair.save_outputs(csv_path, log, list(processed), result, list(paths.values()), before)
            self.assertFalse(csv_path.exists())
            self.assertFalse(log.exists())
            self.assertEqual(repair.fingerprints(list(paths.values())), before)

    def test_schema_errors_are_rejected(self):
        with dataset_fixture() as (root, paths):
            bad = root / "bad.csv"
            for text in ("value,value\n1,2\n", "value\n1,2\n", "other\n1\n", "value\n"):
                with self.subTest(text=text), self.assertRaises(ValueError):
                    bad.write_text(text, encoding="utf-8")
                    repair.load_table(bad, {"value"})


if __name__ == "__main__":
    unittest.main()
