"""Tests for exact-duplicate removal. All file fixtures live outside the repository and are removed.

From the group repository root:
    py -B -m unittest discover -s tests -p "test_remove_duplicates.py" -v
"""

from __future__ import annotations

import copy
import csv
import importlib.util
import io
import json
import runpy
import shutil
import sys
import tempfile
import unittest
import uuid
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from pathlib import Path
from unittest.mock import patch


REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "python scripts/04_remove_duplicates.py"
SPEC = importlib.util.spec_from_file_location("remove_duplicates", SCRIPT)
step4 = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = step4
SPEC.loader.exec_module(step4)
REPAIR = runpy.run_path(str(SCRIPT.with_name("03_repair_solvent_pairs.py")))
COMPOUNDS = [
    {"compound_id": "901", "smiles": "CC(=O)O"},
    {"compound_id": "29", "smiles": "O"},
    {"compound_id": "140", "smiles": "OC(C)C"},
]


def fixture(label="source-Z", value="-0.942900", temperature="283.15", first=0, second=100):
    raw = {"Label": label, "LogS": value, "Temperature, K": temperature,
           "Drug": "Synthetic drug", "SMILES_0": "CC(=O)O", "SMILES_1": "O",
           "SMILES_2": "OC(C)C", "comp_0": "100", "comp_1": str(first), "comp_2": str(second)}
    row = dict(zip(step4.COLUMNS, [value, temperature, "Log solubility", "g/L",
                                 json.dumps([cid for cid, percent in [(29, first), (140, second)] if percent != 0]),
                                 json.dumps([p / 100 for p in [first, second] if p != 0]), "[901]", "[1.0]"]))
    return raw, row


def run_core(records):
    return step4.deduplicate([r for source, r in records], [source["Label"] for source, r in records])


def write_table(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        step4.write_csv(handle, fields, rows)


@contextmanager
def synthetic_dataset(records=None):
    # mkdir inherits Windows temp-directory permissions; unlike Windows 3.13+
    # mkdtemp owner-only ACLs, this also works with a restricted execution token.
    temp_root = Path(tempfile.gettempdir()).resolve()
    folder = temp_root / f"che1147-step4-{uuid.uuid4().hex}"
    folder.mkdir()
    try:
        records = records or [fixture("source-Z"), fixture("source-A", value="0.20"), fixture("source-B")]
        paths = step4.Inputs(folder / "source", folder / step4.INPUT_NAME,
                             folder / "step3_summary.txt", folder / "step3_log.csv",
                             SCRIPT.with_name("02_validate_values.py"), SCRIPT.with_name("03_repair_solvent_pairs.py"))
        raw, rows = [r for r, p in records], [p for r, p in records]
        write_table(paths.raw, list(raw[0]), raw)
        write_table(paths.input_csv, step4.COLUMNS, rows)
        write_table(paths.compounds, ["compound_id", "smiles"], COMPOUNDS)
        paths.metadata.write_text('{"description":"Synthetic source; unit semantics unresolved"}', encoding="utf-8")
        original = []
        for source, row in records:
            former = row.copy()
            former["cmp_ids_solvent"] = "[29, 140]"
            original.append(former)
        changes = REPAIR["repair_rows"](raw, original, COMPOUNDS).changes
        hashes = {step4.relative(p): step4.sha256(p) for p in [paths.raw, paths.compounds, paths.metadata,
                                                             paths.validator, paths.repair_script]}
        entries = [{**c, "input_sha256": json.dumps(hashes)} for c in changes]
        write_table(paths.step3_log, REPAIR["LOG_FIELDS"], entries)
        refresh_evidence(paths)
        outputs = step4.Outputs(folder / "outputs" / step4.DATASET_NAME, folder / "outputs/mapping.csv",
                                folder / "outputs/removals.csv", folder / "outputs/validation.txt")
        yield folder, paths, outputs, records
    finally:
        if folder.resolve().parent != temp_root or not folder.name.startswith("che1147-step4-"):
            raise ValueError("Refusing cleanup outside the synthetic fixture directory")
        shutil.rmtree(folder)


def refresh_evidence(paths):
    paths.step3_summary.write_text("SYNTHETIC STEP 3 EVIDENCE\n\nOUTPUT SHA256\n" +
                                  "\n".join(f"{step4.relative(p)}: {step4.sha256(p)}"
                                            for p in [paths.input_csv, paths.step3_log]) + "\n\n",
                                  encoding="utf-8")


def arguments(paths, records, outputs=None, explicit_dry=False):
    result = run_core(records)
    args = ["--input-csv", str(paths.input_csv), "--dataset-dir", str(paths.dataset_dir),
            "--step3-summary", str(paths.step3_summary), "--step3-log", str(paths.step3_log),
            "--expected-input-records", str(len(records)),
            "--expected-duplicate-groups", str(len({e['duplicate_group'] for e in result.log})),
            "--expected-removed-records", str(len(result.log)), "--expected-output-records", str(len(result.rows)),
            "--expected-retained-labels", *[e["retained_raw_label"] for e in result.log],
            "--expected-removed-labels", *[e["removed_raw_label"] for e in result.log]]
    if outputs is not None:
        for flag, path in zip(["--output-csv", "--row-mapping", "--removal-log", "--validation-summary"], outputs.paths()):
            args += [flag, str(path)]
    if explicit_dry:
        args.append("--dry-run")
    return args


def call_main(args):
    out, err = io.StringIO(), io.StringIO()
    with redirect_stdout(out), redirect_stderr(err):
        code = step4.main(args)
    return code, out.getvalue(), err.getvalue()


class ExactRowPolicyTests(unittest.TestCase):
    def test_no_duplicates(self):
        records = [fixture("z", "-3"), fixture("a", "2")]
        result = run_core(records)
        self.assertEqual(result.rows, [row for raw, row in records])
        self.assertEqual(result.log, [])
        self.assertEqual([m["new_data_row"] for m in result.mapping], ["1", "2"])

    def test_pair_keeps_first_and_preserves_original_label(self):
        records = [fixture("z"), fixture("a")]
        result = run_core(records)
        self.assertEqual(len(result.rows), 1)
        self.assertEqual(result.log[0]["removed_raw_label"], "a")
        self.assertEqual(result.log[0]["retained_raw_label"], "z")
        self.assertEqual(result.log[0]["reason"], step4.REASON)

    def test_more_than_two_occurrences_share_one_group(self):
        result = run_core([fixture("first"), fixture("second"), fixture("third"), fixture("fourth")])
        self.assertEqual(len(result.rows), 1)
        self.assertEqual(len(result.log), 3)
        self.assertEqual({e["duplicate_group"] for e in result.log}, {"exact_duplicate_001"})
        self.assertEqual({e["retained_original_index"] for e in result.log}, {"0"})

    def test_separated_duplicates_shift_positions_without_sorting(self):
        records = [fixture("z", "2"), fixture("y", "1"), fixture("x", "2"), fixture("w", "3")]
        result = run_core(records)
        self.assertEqual([r["value"] for r in result.rows], ["2", "1", "3"])
        self.assertEqual(result.mapping[2], dict(zip(step4.MAPPING_FIELDS, ["3", "", "x", "removed_duplicate", "2", ""])))
        self.assertEqual(result.mapping[3]["new_index"], "2")
        self.assertEqual(result.log[0]["retained_new_data_row"], "1")
        step4.validate_mapping([p for r, p in records], [r for r, p in records], result)

    def test_identical_conditions_different_solubility_remain(self):
        records = [fixture("a", "1.00000000000000001"), fixture("b", "1.00000000000000002")]
        result = run_core(records)
        self.assertEqual(len(result.rows), 2)
        stats = step4.output_statistics(result.rows)
        self.assertEqual(stats["repeated_condition_groups"], 1)
        self.assertEqual(stats["different_solubility_condition_groups"], 1)

    def test_no_numeric_or_list_text_normalization_in_duplicate_key(self):
        records = [fixture("a", "0.10"), fixture("b", "0.1"), fixture("c", "0.10")]
        records[2][1]["cmp_ids_drug"] = "[901.0]"
        result = run_core(records)
        self.assertEqual(len(result.rows), 3)
        self.assertEqual(result.log, [])

    def test_every_column_participates_without_component_reordering(self):
        a, b = fixture("a", first=25, second=75), fixture("b", first=25, second=75)
        b[1]["cmp_ids_solvent"] = "[140, 29]"
        b[1]["cmp_mole_fractions_solvent"] = "[0.75, 0.25]"
        self.assertEqual(len(run_core([a, b]).rows), 2)
        for column in step4.COLUMNS:
            with self.subTest(column=column):
                records = [fixture("a"), fixture("b")]
                records[1][1][column] += " "
                self.assertEqual(len(run_core(records).rows), 2)

    def test_no_input_mutation_and_all_cell_strings_preserved(self):
        records = [fixture("z"), fixture("a"), fixture("b", "1.0")]
        records[0][1]["unit"] = records[1][1]["unit"] = 'custom, "quoted"\nunit '
        records[0][1]["cmp_mole_fractions_solvent"] = records[1][1]["cmp_mole_fractions_solvent"] = "[1.000000]"
        before = copy.deepcopy(records)
        result = run_core(records)
        self.assertEqual(records, before)
        self.assertEqual(result.rows, [before[0][1], before[2][1]])
        self.assertEqual(list(result.rows[0]), step4.COLUMNS)

    def test_second_pass_removes_nothing(self):
        records = [fixture("a"), fixture("b"), fixture("c", "2")]
        first = run_core(records)
        labels = [m["raw_label"] for m in first.mapping if m["status"] == "retained"]
        second = step4.deduplicate(first.rows, labels)
        self.assertEqual(second.rows, first.rows)
        self.assertEqual(second.log, [])

    def test_rejects_unverified_or_ambiguous_labels_and_invalid_row_schema(self):
        row = fixture()[1]
        for rows, labels in [([row], []), ([row, row], ["same", "same"]), ([row], [""]),
                             ([{**row, "index": "0"}], ["source"])]:
            with self.subTest(labels=labels), self.assertRaises(ValueError):
                step4.deduplicate(rows, labels)

    def test_mapping_validation_rejects_cell_link_position_and_label_tampering(self):
        records = [fixture("a"), fixture("b"), fixture("c", "2")]
        for mutate in [lambda r: r.rows[0].update(value="99"),
                       lambda r: r.mapping[1].update(new_data_row="0"),
                       lambda r: r.mapping[2].update(new_index="4"),
                       lambda r: r.mapping[0].update(raw_label="invented"),
                       lambda r: r.log[0].update(retained_original_index="2")]:
            result = run_core(records)
            mutate(result)
            with self.assertRaises(ValueError):
                step4.validate_mapping([p for r, p in records], [r for r, p in records], result)


class FilesAndValidationTests(unittest.TestCase):
    def test_default_and_explicit_dry_runs_write_no_files_or_directories(self):
        with synthetic_dataset() as (folder, paths, outputs, records):
            before = {p.relative_to(folder): step4.sha256(p) for p in folder.rglob("*") if p.is_file()}
            for explicit in [False, True]:
                code, out, err = call_main(arguments(paths, records, explicit_dry=explicit))
                self.assertEqual(code, 0, err)
                self.assertIn("DRY RUN", out)
                self.assertFalse(outputs.dataset.parent.exists())
            after = {p.relative_to(folder): step4.sha256(p) for p in folder.rglob("*") if p.is_file()}
            self.assertEqual(before, after)

    def test_saved_synthetic_outputs_are_read_back_and_quoted_cells_preserved(self):
        records = [fixture("z", first=25, second=75), fixture("x", "1"), fixture("a", first=25, second=75)]
        records[0][1]["unit"] = records[2][1]["unit"] = 'g/L, "quoted"\nlabel'
        with synthetic_dataset(records) as (folder, paths, outputs, records):
            reads = []
            real_read = step4.read_csv

            def trace(path, **kwargs):
                reads.append(path.resolve())
                return real_read(path, **kwargs)

            with patch.object(step4, "read_csv", side_effect=trace):
                code, out, err = call_main(arguments(paths, records, outputs))
            self.assertEqual(code, 0, err)
            self.assertIn("SAVED OUTPUT VALIDATION", out)
            for path in outputs.paths()[:3]:
                self.assertIn(path.resolve(), reads)
            rows = real_read(outputs.dataset, fields=step4.COLUMNS)
            self.assertEqual(rows, [records[0][1], records[1][1]])
            mapping = real_read(outputs.mapping, fields=step4.MAPPING_FIELDS)
            self.assertEqual(mapping[2]["new_data_row"], "")
            self.assertEqual(mapping[2]["new_index"], "")
            log = real_read(outputs.log, fields=step4.LOG_FIELDS)
            self.assertEqual(log[0]["retained_raw_label"], "z")
            summary = outputs.summary.read_text(encoding="utf-8")
            self.assertIn(step4.DATASET_NAME, summary)
            self.assertIn("positional indices cannot be reused unchanged", summary)
            self.assertNotIn(str(REPO_ROOT), summary)
            self.assertNotIn(str(folder), summary)
            for path in outputs.paths()[:3]:
                self.assertIn(step4.sha256(path), summary)

    def test_invalid_schemas_and_csv_width_or_quoting_are_rejected(self):
        with synthetic_dataset() as (folder, paths, outputs, records):
            for text in ["value,value\n1,2\n", ",".join(reversed(step4.COLUMNS)) + "\n",
                         ",".join(step4.COLUMNS) + "\n1,2\n",
                         ",".join(step4.COLUMNS) + '\n"unclosed']:
                paths.input_csv.write_text(text, encoding="utf-8")
                with self.assertRaises((ValueError, csv.Error)):
                    step4.read_csv(paths.input_csv, fields=step4.COLUMNS)

    def test_step3_hash_mismatch_stops_before_saving(self):
        with synthetic_dataset() as (folder, paths, outputs, records):
            paths.input_csv.write_bytes(paths.input_csv.read_bytes() + b"\n")
            code, out, err = call_main(arguments(paths, records, outputs))
            self.assertEqual(code, 2)
            self.assertFalse(outputs.dataset.parent.exists())

    def test_unverified_target_temperature_drug_and_fraction_correspondence_stops(self):
        for field, value in [("value", "999"), ("Temperature, K", "1"),
                             ("cmp_ids_drug", "[29]"), ("cmp_ids_solvent", "[29]")]:
            with self.subTest(field=field), synthetic_dataset() as (folder, paths, outputs, records):
                rows = [p.copy() for r, p in records]
                rows[0][field] = value
                write_table(paths.input_csv, step4.COLUMNS, rows)
                refresh_evidence(paths)
                code, out, err = call_main(arguments(paths, records, outputs))
                self.assertEqual(code, 2)
                self.assertFalse(outputs.dataset.parent.exists())

    def test_row_misalignment_stops_even_with_updated_output_hash(self):
        with synthetic_dataset() as (folder, paths, outputs, records):
            write_table(paths.input_csv, step4.COLUMNS, [records[1][1], records[0][1], records[2][1]])
            refresh_evidence(paths)
            self.assertEqual(call_main(arguments(paths, records, outputs))[0], 2)
            self.assertFalse(outputs.dataset.parent.exists())

    def test_ambiguous_compound_reference_is_rejected(self):
        with synthetic_dataset() as (folder, paths, outputs, records):
            write_table(paths.compounds, ["compound_id", "smiles"], COMPOUNDS + [COMPOUNDS[0]])
            self.assertEqual(call_main(arguments(paths, records, outputs))[0], 2)
            self.assertFalse(outputs.dataset.parent.exists())

    def test_step3_log_field_and_protected_fingerprint_tampering_are_rejected(self):
        for field, value in [("raw_label", "wrong"), ("comp_2", "25"), ("input_sha256", '{}')]:
            with self.subTest(field=field), synthetic_dataset() as (folder, paths, outputs, records):
                log = step4.read_csv(paths.step3_log, fields=REPAIR["LOG_FIELDS"])
                log[0][field] = value
                write_table(paths.step3_log, REPAIR["LOG_FIELDS"], log)
                refresh_evidence(paths)
                self.assertEqual(call_main(arguments(paths, records, outputs))[0], 2)
                self.assertFalse(outputs.dataset.parent.exists())

    def test_original_processed_csv_is_not_a_cleaning_input(self):
        with synthetic_dataset() as (folder, paths, outputs, records):
            forbidden = paths.dataset_dir / "processed_data/processed_DrugSolubilityData.csv"
            write_table(forbidden, step4.COLUMNS, [p for r, p in records])
            paths.input_csv = forbidden
            refresh_evidence(paths)
            self.assertEqual(call_main(arguments(paths, records, outputs))[0], 2)
            self.assertFalse(outputs.dataset.parent.exists())

    def test_input_and_existing_output_collisions_do_not_change_existing_files(self):
        with synthetic_dataset() as (folder, paths, outputs, records):
            before = step4.sha256(paths.input_csv)
            blocked = copy.copy(outputs)
            blocked.dataset = paths.input_csv
            self.assertEqual(call_main(arguments(paths, records, blocked))[0], 2)
            self.assertEqual(step4.sha256(paths.input_csv), before)
            for index, path in enumerate(outputs.paths()):
                with self.subTest(index=index):
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_text("existing sentinel", encoding="utf-8")
                    self.assertEqual(call_main(arguments(paths, records, outputs))[0], 2)
                    self.assertEqual(path.read_text(encoding="utf-8"), "existing sentinel")
                    path.unlink()

    def test_output_destinations_must_be_distinct(self):
        with synthetic_dataset() as (folder, paths, outputs, records):
            outputs.log = outputs.mapping
            self.assertEqual(call_main(arguments(paths, records, outputs))[0], 2)
            self.assertFalse(outputs.dataset.parent.exists())

    def test_all_four_save_arguments_required_and_explicit_dry_conflicts(self):
        with synthetic_dataset() as (folder, paths, outputs, records):
            for args in [arguments(paths, records) + ["--output-csv", str(outputs.dataset)],
                         arguments(paths, records, outputs, explicit_dry=True)]:
                with self.assertRaises(SystemExit) as raised, redirect_stderr(io.StringIO()):
                    step4.main(args)
                self.assertEqual(raised.exception.code, 2)
                self.assertFalse(outputs.dataset.parent.exists())

    def test_unexpected_counts_or_label_identities_prevent_saving(self):
        with synthetic_dataset() as (folder, paths, outputs, records):
            base = arguments(paths, records, outputs)
            for extra in [["--expected-removed-records", "99"], ["--expected-retained-labels", "invented"]]:
                self.assertEqual(call_main(base + extra)[0], 2)
                self.assertFalse(outputs.dataset.parent.exists())

    def test_saved_validation_failure_is_reported_without_claiming_success(self):
        with synthetic_dataset() as (folder, paths, outputs, records):
            real_validate = step4.validate_result

            def fail_saved(*args, **kwargs):
                if kwargs["saved"]:
                    raise ValueError("synthetic readback failure")
                return real_validate(*args, **kwargs)

            with patch.object(step4, "validate_result", side_effect=fail_saved):
                code, out, err = call_main(arguments(paths, records, outputs))
            self.assertEqual(code, 2)
            self.assertIn("Created outputs left for review", err)
            self.assertIn("synthetic readback failure", err)
            self.assertNotIn("SAVED OUTPUT VALIDATION", out)
            self.assertEqual(outputs.summary.read_text(encoding="utf-8"), "")

    def test_saved_mapping_validation_uses_raw_mapping_after_shift(self):
        records = [fixture("z", "1"), fixture("y", "1"), fixture("x", "2", temperature="300.15")]
        with synthetic_dataset(records) as (folder, paths, outputs, records):
            code, out, err = call_main(arguments(paths, records, outputs))
            self.assertEqual(code, 0, err)
            self.assertEqual(step4.read_csv(outputs.dataset, fields=step4.COLUMNS)[1]["value"], "2")
            mappings = step4.read_csv(outputs.mapping, fields=step4.MAPPING_FIELDS)
            self.assertEqual(mappings[2]["raw_label"], "x")
            self.assertEqual(mappings[2]["new_data_row"], "2")


if __name__ == "__main__":
    unittest.main()
