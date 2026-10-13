#!/usr/bin/env python3
"""Repair only solvent IDs left behind when zero fractions were dropped.

Python 3.10+; standard library only. From the group repository root:
    py -B "python scripts/03_repair_solvent_pairs.py" --dry-run

The default dataset is the sibling chemixhub/datasets/drug-solubility folder,
resolved from this script's location, independently of the working directory.
No output is written unless BOTH --output-csv and --change-log are provided.
Those options cannot be combined with --dry-run. Existing outputs and inputs
are never overwritten. --processed-input can select a previously repaired CSV.

Only cmp_ids_solvent may change. All other cells, row order, duplicates, unit
labels, and record count are preserved. Fractions are not rounded or normalized.
Raw SMILES are mapped using explicit compound_id values from compounds.csv.
Records are compared by position after checking their contents; Label is only
used in diagnostics. Unexplained differences abort the entire operation.

Exit codes: 0 = original input checked and corrected rows pass validation;
            2 = input, correspondence, validation, or output-path failure.
Review findings from Step 2 are reported, preserved, and not resolved here.
"""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import math
import runpy
import sys
from collections import Counter
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path
from types import FunctionType


REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DATASET = REPO_ROOT.parent / "chemixhub" / "datasets" / "drug-solubility"
VALIDATOR_PATH = Path(__file__).with_name("02_validate_values.py")
SUM_TOL = 1e-6
PAIR_TOL = 1e-12  # Floating-point serialization only; never used to detect zero.
ID_COLUMN = "cmp_ids_solvent"
RAW_REQUIRED = {
    "LogS", "Temperature, K", "Drug", "SMILES_0", "SMILES_1", "SMILES_2",
    "comp_0", "comp_1", "comp_2",
}
PROCESSED_REQUIRED = {
    "value", "Temperature, K", "property", "unit", ID_COLUMN,
    "cmp_mole_fractions_solvent", "cmp_ids_drug", "cmp_mole_fractions_drug",
}
LOG_FIELDS = [
    "data_row", "raw_label", "zero_solvent_position", "changed_column",
    "before", "after", "cmp_mole_fractions_solvent", "SMILES_1", "comp_1",
    "SMILES_2", "comp_2", "input_sha256",
]


@dataclass
class Table:
    header: list[str]
    rows: list[dict[str, str]]


@dataclass
class RepairResult:
    rows: list[dict[str, str]]
    changes: list[dict]


def load_table(path: Path, required: set[str]) -> Table:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        header = reader.fieldnames or []
        if not header or len(header) != len(set(header)):
            raise ValueError(f"{path.name}: empty or duplicate column names")
        if required - set(header):
            raise ValueError(f"{path.name}: missing columns {sorted(required - set(header))}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{path.name}: no data records")
    if any(None in row or None in row.values() for row in rows):
        raise ValueError(f"{path.name}: incorrect record width")
    return Table(header, rows)


def finite_number(value: object, context: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{context}: boolean is not a number")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{context}: invalid number {value!r}") from exc
    if not math.isfinite(result):
        raise ValueError(f"{context}: nonfinite number")
    return result


def compound_mapping(rows: list[dict[str, str]]) -> dict[str, int]:
    mapping = {}
    known_ids = set()
    for row in rows:
        cid = finite_number(row["compound_id"], "compounds.csv compound_id")
        smiles = row["smiles"].strip()
        if cid < 0 or not cid.is_integer() or not smiles:
            raise ValueError("compounds.csv: invalid ID or empty SMILES")
        cid = int(cid)
        if cid in known_ids or smiles in mapping:
            raise ValueError("compounds.csv: duplicate ID or ambiguous SMILES mapping")
        mapping[smiles] = cid
        known_ids.add(cid)
    if not mapping:
        raise ValueError("compounds.csv: empty mapping")
    return mapping


def numeric_list(text: str, context: str, known_ids: set[int] | None = None) -> list:
    try:
        if len(text) > 10000:
            raise ValueError("oversized list")
        items = ast.literal_eval(text)
    except (ValueError, SyntaxError, TypeError, RecursionError) as exc:
        raise ValueError(f"{context}: invalid list") from exc
    if not isinstance(items, list) or not items:
        raise ValueError(f"{context}: expected a nonempty list")
    for item in items:
        if isinstance(item, bool) or not isinstance(item, (int, float)):
            raise ValueError(f"{context}: invalid numeric list item")
        finite_number(item, context)
    if known_ids is not None:
        if any(x < 0 or x != int(x) or int(x) not in known_ids for x in items):
            raise ValueError(f"{context}: invalid or unknown compound ID")
        if len(set(items)) != len(items):
            raise ValueError(f"{context}: repeated compound ID")
        return [int(x) for x in items]
    if any(x < 0 or x > 1 for x in items) or not math.isclose(
        sum(items), 1, rel_tol=0, abs_tol=SUM_TOL
    ):
        raise ValueError(f"{context}: invalid fractions or fraction sum")
    return items


def same_fractions(left: list, right: list) -> bool:
    return len(left) == len(right) and all(
        math.isclose(x, y, rel_tol=0, abs_tol=PAIR_TOL)
        for x, y in zip(left, right)
    )


def repair_rows(raw: list[dict[str, str]], processed: list[dict[str, str]],
                compounds: list[dict[str, str]]) -> RepairResult:
    """Check every row and return copies; never mutate any supplied record."""
    if len(raw) != len(processed):
        raise ValueError("Raw/processed row count mismatch; cannot establish correspondence")
    mapping = compound_mapping(compounds)
    known_ids = set(mapping.values())
    corrected = []
    changes = []
    for row_n, (source, original) in enumerate(zip(raw, processed), start=1):
        context = f"data row {row_n} (raw Label {source.get('Label', '')})"
        for raw_key, processed_key in [("LogS", "value"),
                                       ("Temperature, K", "Temperature, K")]:
            left = finite_number(source[raw_key], f"{context} raw {raw_key}")
            right = finite_number(original[processed_key], f"{context} {processed_key}")
            if left != right or (raw_key == "Temperature, K" and left <= 0):
                raise ValueError(f"{context}: correspondence failure for {processed_key}")
        if original["property"].strip() != "Log solubility" or not original["unit"].strip():
            raise ValueError(f"{context}: unexpected property or empty unit")
        percents = [finite_number(source[f"comp_{k}"], f"{context} comp_{k}")
                    for k in range(3)]
        if (any(x < 0 or x > 100 for x in percents)
                or not math.isclose(percents[0], 100, rel_tol=0, abs_tol=100 * SUM_TOL)
                or not math.isclose(sum(percents[1:]), 100, rel_tol=0,
                                    abs_tol=100 * SUM_TOL)):
            raise ValueError(f"{context}: invalid raw composition; no normalization performed")
        raw_ids = [mapping.get(source[f"SMILES_{k}"].strip()) for k in range(3)]
        if any(cid is None and percent > 0 for cid, percent in zip(raw_ids, percents)):
            raise ValueError(f"{context}: positive raw component has no SMILES mapping")
        drug_ids = numeric_list(original["cmp_ids_drug"], f"{context} drug IDs", known_ids)
        drug_fractions = numeric_list(original["cmp_mole_fractions_drug"],
                                      f"{context} drug fractions")
        if drug_ids != [raw_ids[0]] or not same_fractions(drug_fractions, [percents[0] / 100]):
            raise ValueError(f"{context}: drug correspondence failure")

        # Pair using raw positions FIRST. Only an exactly-zero component is absent.
        positive_pairs = [(raw_ids[k], percents[k] / 100)
                          for k in (1, 2) if percents[k] > 0]
        expected_ids = [cid for cid, fraction in positive_pairs]
        expected_fractions = [fraction for cid, fraction in positive_pairs]
        ids = numeric_list(original[ID_COLUMN], f"{context} solvent IDs", known_ids)
        fractions = numeric_list(original["cmp_mole_fractions_solvent"],
                                 f"{context} solvent fractions")
        output = original.copy()
        if len(ids) == len(fractions):
            # Preserve already-correct records, including explicit zero pairs.
            positive_match = ids == expected_ids and same_fractions(fractions, expected_fractions)
            full_match = ids == raw_ids[1:] and same_fractions(
                fractions, [x / 100 for x in percents[1:]]
            )
            if not (positive_match or full_match):
                raise ValueError(f"{context}: solvent pairs do not match raw positions/fractions")
        else:
            zero_positions = [k for k in (1, 2) if percents[k] == 0]
            if len(ids) != 2 or len(fractions) != 1 or len(zero_positions) != 1:
                raise ValueError(f"{context}: length mismatch does not fit the zero-fraction defect")
            if None in raw_ids[1:] or ids != raw_ids[1:]:
                raise ValueError(f"{context}: solvent ID correspondence failure; will not guess")
            if not same_fractions(fractions, expected_fractions):
                raise ValueError(f"{context}: existing fraction list differs from raw nonzero fractions")
            output[ID_COLUMN] = json.dumps(expected_ids)
            changes.append({
                "data_row": row_n, "raw_label": source.get("Label", ""),
                "zero_solvent_position": zero_positions[0], "changed_column": ID_COLUMN,
                "before": original[ID_COLUMN], "after": output[ID_COLUMN],
                "cmp_mole_fractions_solvent": original["cmp_mole_fractions_solvent"],
                **{key: source[key] for key in ("SMILES_1", "comp_1", "SMILES_2", "comp_2")},
            })
        corrected.append(output)
    result = RepairResult(corrected, changes)
    check_preservation(processed, result)
    return result


def check_preservation(original: list[dict[str, str]], result: RepairResult) -> None:
    if len(original) != len(result.rows):
        raise ValueError("Internal preservation failure: record count changed")
    changed_rows = {change["data_row"] for change in result.changes}
    for row_n, (before, after) in enumerate(zip(original, result.rows), start=1):
        differing = {key for key in before.keys() | after.keys() if before.get(key) != after.get(key)}
        if differing != ({ID_COLUMN} if row_n in changed_rows else set()):
            raise ValueError(f"Internal preservation failure at data row {row_n}")


def input_paths(dataset_dir: Path, processed_input: Path | None = None) -> dict[str, Path]:
    original_processed = dataset_dir / "processed_data" / "processed_DrugSolubilityData.csv"
    return {
        "raw": dataset_dir / "raw_data" / "DrugSolubilityData.csv",
        "processed": processed_input if processed_input is not None else original_processed,
        "original_processed": original_processed,
        "compounds": dataset_dir / "processed_data" / "compounds.csv",
        "metadata": dataset_dir / "processed_data" / "croissant.json",
        "upstream_processor": dataset_dir / "processed_data" / "data_processing.py",
        "validator": VALIDATOR_PATH,
    }


def fingerprints(paths: list[Path]) -> dict[str, str]:
    return {str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in paths if path.is_file()}


def load_validator(path: Path) -> dict:
    # run_path compiles source directly, without editing it or creating a .pyc.
    return runpy.run_path(str(path))


def validate_with_existing(namespace: dict, paths: dict[str, Path],
                           memory_rows: list[dict[str, str]] | None = None) -> tuple:
    """Use unchanged Step 2 checks, substituting only the processed rows in memory.

    Copy the function globals rather than patching the validator's namespace.
    Its summary hashes always describe disk inputs, even for a memory check;
    the CLI therefore prints distinct counts instead of its full summary.
    """
    original = namespace["validate"]
    check = original
    if memory_rows is not None:
        environment = original.__globals__.copy()
        disk_loader = environment["load_csv"]

        def memory_loader(path: Path, required: set[str]) -> list[dict[str, str]]:
            if path.resolve() == paths["processed"].resolve():
                if any(required - set(row) for row in memory_rows):
                    raise ValueError("Corrected records are missing required columns")
                return [row.copy() for row in memory_rows]
            return disk_loader(path, required)

        environment["load_csv"] = memory_loader
        check = FunctionType(original.__code__, environment, original.__name__,
                             original.__defaults__, original.__closure__)
    return check(paths["raw"], paths["processed"], paths["compounds"], paths["metadata"], SUM_TOL)


def print_validation(label: str, issues: list[dict]) -> None:
    errors = Counter(issue["check"] for issue in issues if issue["severity"] == "ERROR")
    reviews = [issue for issue in issues if issue["severity"] == "REVIEW"]
    print(f"\n{label}")
    print(f"Structural errors: {sum(errors.values()):,}")
    for check, count in sorted(errors.items()):
        print(f"  {check}: {count:,}")
    print(f"Review items: {len(reviews)}")
    for issue in reviews:
        print(f"  {issue['check']}: {issue['detail']}")


def guard_outputs(output_csv: Path, change_log: Path, protected: list[Path]) -> None:
    destinations = [output_csv.resolve(), change_log.resolve()]
    if destinations[0] == destinations[1]:
        raise ValueError("Repaired CSV and change log must have different paths")
    protected_paths = {path.resolve() for path in protected}
    for supplied, resolved in zip((output_csv, change_log), destinations):
        if resolved in protected_paths:
            raise ValueError(f"Output would overwrite a protected input: {supplied}")
        if supplied.exists() or supplied.is_symlink():
            raise ValueError(f"Output already exists; refusing to overwrite: {supplied}")


def save_outputs(output_csv: Path, change_log: Path, header: list[str], result: RepairResult,
                 protected: list[Path], hashes: dict[str, str]) -> None:
    """Create both files exclusively, after all checks; remove partial new files on failure."""
    guard_outputs(output_csv, change_log, protected)
    created = []
    try:
        with ExitStack() as stack:
            handles = []
            for path in (output_csv, change_log):
                path.parent.mkdir(parents=True, exist_ok=True)
                handle = stack.enter_context(path.open("x", encoding="utf-8", newline=""))
                created.append(path)
                handles.append(handle)
            writer = csv.DictWriter(handles[0], fieldnames=header)
            writer.writeheader()
            writer.writerows(result.rows)
            log_writer = csv.DictWriter(handles[1], fieldnames=LOG_FIELDS)
            log_writer.writeheader()
            for change in result.changes:
                log_writer.writerow({**change, "input_sha256": json.dumps(hashes, sort_keys=True)})
    except BaseException:
        for path in created:
            path.unlink(missing_ok=True)
        raise


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dataset-dir", type=Path, default=DEFAULT_DATASET)
    parser.add_argument("--processed-input", type=Path,
                        help="Optional processed CSV, e.g. a previously repaired copy")
    parser.add_argument("--dry-run", action="store_true", help="Console only; this is also the default")
    parser.add_argument("--output-csv", type=Path, help="Explicit NEW repaired CSV destination")
    parser.add_argument("--change-log", type=Path, help="Explicit NEW change-log CSV destination")
    args = parser.parse_args(argv)
    if (args.output_csv is None) != (args.change_log is None):
        parser.error("--output-csv and --change-log must be supplied together")
    if args.dry_run and args.output_csv is not None:
        parser.error("--dry-run cannot be combined with output arguments")
    paths = input_paths(args.dataset_dir, args.processed_input)
    protected = list(paths.values()) + [Path(__file__)]
    try:
        if args.output_csv is not None:
            guard_outputs(args.output_csv, args.change_log, protected)
        hashes = fingerprints(protected)
        raw = load_table(paths["raw"], RAW_REQUIRED)
        processed = load_table(paths["processed"], PROCESSED_REQUIRED)
        compounds = load_table(paths["compounds"], {"compound_id", "smiles"})
        namespace = load_validator(paths["validator"])
        _, original_issues, original_candidates = validate_with_existing(namespace, paths)
        print("STEP 3 - SOLVENT ID REPAIR")
        print("Mode: DRY RUN (no output files)" if args.output_csv is None else "Mode: explicit separate outputs")
        print(f"Raw records: {len(raw.rows):,}; processed records: {len(processed.rows):,}")
        print_validation("ORIGINAL DATA VALIDATION (unchanged disk input)", original_issues)
        unexpected = [issue for issue in original_issues
                      if issue["severity"] == "ERROR" and issue["check"] != "SOLVENT_LENGTH_MISMATCH"]
        if unexpected:
            first = unexpected[0]
            raise ValueError(f"Original data has an unrelated error: {first['check']} "
                             f"at data row {first['data_row']}: {first['detail']}")
        result = repair_rows(raw.rows, processed.rows, compounds.rows)
        original_errors = sum(issue["severity"] == "ERROR" for issue in original_issues)
        if len(result.changes) != original_errors:
            raise ValueError("Repair count does not agree with original validation findings")
        _, corrected_issues, corrected_candidates = validate_with_existing(namespace, paths, result.rows)
        print_validation("CORRECTED DATA VALIDATION (in-memory records; no repaired file read)", corrected_issues)
        if any(issue["severity"] == "ERROR" for issue in corrected_issues):
            raise ValueError("Corrected records failed validation; no output will be saved")
        if (original_candidates != corrected_candidates
                or [i for i in original_issues if i["severity"] == "REVIEW"]
                != [i for i in corrected_issues if i["severity"] == "REVIEW"]):
            raise ValueError("Review findings changed unexpectedly")
        second = repair_rows(raw.rows, result.rows, compounds.rows)
        if second.changes or second.rows != result.rows:
            raise ValueError("Second repair pass changed the corrected records")
        if fingerprints(protected) != hashes:
            raise ValueError("An input or script changed during processing; output refused")
        zero_counts = Counter(change["zero_solvent_position"] for change in result.changes)
        print("\nREPAIR SUMMARY")
        print(f"Affected records detected from contents: {len(result.changes):,}")
        print(f"Solvent 1 zero: {zero_counts[1]:,}; solvent 2 zero: {zero_counts[2]:,}")
        print(f"Records retained in original order: {len(result.rows):,}")
        print(f"Unchanged records: {len(result.rows) - len(result.changes):,}")
        print("Only changed column: cmp_ids_solvent (all other cells preserved)")
        print("Second pass changes: 0")
        print(f"Repeated-condition candidate rows retained: {len(corrected_candidates):,}")
        print("Input and script SHA256 check: unchanged")
        for path, digest in hashes.items():
            print(f"  {path}: {digest}")
        if args.output_csv is not None:
            save_outputs(args.output_csv, args.change_log, processed.header, result, protected, hashes)
            print(f"Repaired CSV written to: {args.output_csv.resolve()}")
            print(f"Change log written to: {args.change_log.resolve()}")
        else:
            print("No repaired dataset or change log saved.")
        return 0
    except (OSError, ValueError, csv.Error) as exc:
        print(f"Repair stopped: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
