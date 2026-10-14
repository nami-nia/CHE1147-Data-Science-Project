#!/usr/bin/env python3
"""Select binary-solvent mixtures from the existing three-drug subset.

Research-scope filtering only; no cleaning, sorting, deduplication or outlier
removal. Python 3.10+; standard library only. From the repository root:

    py -B "python scripts/06_remove_pure_solvents.py"
    py -B "python scripts/06_remove_pure_solvents.py" --input "data/processed/subset_drugs_DrugSolubilityData.csv" --output "data/processed/subset_drugs_binary_only_DrugSolubilityData.csv"

Relative arguments start at this script's repository, regardless of the current
working directory. Default discovery requires exactly one existing input named
subset_drugs_DrugSolubilityData.csv. The output defaults beside that input when
inside the repository, otherwise to data/processed/. Outputs stay inside the
repository and use subset_drugs_binary_only_DrugSolubilityData.csv.

Every record must have one drug ID in {5, 43, 86}, distinct solvent IDs, matching
nonempty fraction lists, and finite fractions in [0, 1] summing to 1 with absolute
tolerance 1e-6 and zero relative tolerance. Exactly positive fractions determine
presence: retain two active solvents; exclude one; reject any other composition.
Expected counts only verify the audited result; they never select records.
All decoded CSV fields, columns and retained-record order remain unchanged.

A new output is written temporarily, read back completely and published without
overwriting only after all checks pass. The published file is also read back.
Matching existing outputs are verified and reused unchanged; conflicting outputs
are refused. No overwrite option, reports, permanent test fixtures or Git changes.
Parser/filter self-checks run in memory. Exit codes: 0 = verified; 2 = failure.
"""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import math
import os
import sys
import uuid
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
INPUT_NAME = "subset_drugs_DrugSolubilityData.csv"
OUTPUT_NAME = "subset_drugs_binary_only_DrugSolubilityData.csv"
DRUG_COLUMN = "cmp_ids_drug"
SOLVENT_COLUMN = "cmp_ids_solvent"
FRACTION_COLUMN = "cmp_mole_fractions_solvent"
EXPECTED_COUNTS = {5: (390, 30, 360), 43: (189, 27, 162), 86: (105, 15, 90)}
DRUG_NAMES = {5: "Itraconazole", 43: "Posaconazole", 86: "Ketoconazole"}
SUM_TOLERANCE = 1e-6


@dataclass
class FilterResult:
    retained: list[list[str]]
    retained_positions: list[int]
    excluded_positions: list[int]
    counts: dict[int, dict[str, int]]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def resolve_path(path: Path) -> Path:
    return (path if path.is_absolute() else REPO_ROOT / path).resolve()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def locate_input(requested: Path | None) -> Path:
    if requested is None:
        copies = [p for p in REPO_ROOT.rglob(INPUT_NAME)
                  if p.is_file() and ".git" not in p.relative_to(REPO_ROOT).parts]
        require(bool(copies), f"Required input {INPUT_NAME} is missing in {REPO_ROOT}")
        require(len(copies) == 1, "Multiple input copies found; choose explicitly with --input:\n"
                + "\n".join(str(p) for p in sorted(copies)))
        requested = copies[0]
    path = resolve_path(requested)
    require(path.name == INPUT_NAME, f"Input filename must be {INPUT_NAME}: {path}")
    require(path.is_file(), f"Required input is missing: {path}")
    return path


def validate_header(header: list[str]) -> dict[str, int]:
    require(bool(header) and all(name.strip() for name in header), "Empty or unnamed CSV columns")
    require(len(header) == len(set(header)), "Ambiguous CSV header: duplicate column names")
    for name in (DRUG_COLUMN, SOLVENT_COLUMN, FRACTION_COLUMN):
        require(header.count(name) == 1, f"CSV header must contain {name} exactly once")
    return {name: header.index(name) for name in (DRUG_COLUMN, SOLVENT_COLUMN, FRACTION_COLUMN)}


def read_csv(path: Path) -> tuple[list[str], list[list[str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, strict=True)
        header = next(reader, [])
        validate_header(header)
        rows = []
        try:
            for number, row in enumerate(reader, 1):
                require(len(row) == len(header),
                        f"{path}: data record {number} has {len(row)} fields; expected {len(header)}")
                rows.append(row)
        except csv.Error as exc:
            raise ValueError(f"{path}: invalid CSV quoting at data record {len(rows) + 1}: {exc}") from exc
    return header, rows


def parse_list(cell: str, column: str, number: int) -> list:
    try:
        parsed = ast.literal_eval(cell)
        require(type(parsed) is list and bool(parsed), "expected a nonempty list")
        return parsed
    except (ValueError, SyntaxError, TypeError, RecursionError, MemoryError, OverflowError) as exc:
        raise ValueError(f"Data record {number}: invalid {column} cell {cell!r}: {exc}") from exc


def parse_ids(cell: str, column: str, number: int) -> list[int]:
    values = parse_list(cell, column, number)
    for value in values:
        require(type(value) in (int, float) and value >= 0
                and (type(value) is int or (math.isfinite(value) and value.is_integer())),
                f"Data record {number}: invalid {column} cell {cell!r}; "
                "IDs must be non-negative integers or finite, exactly integral floats")
    return [int(value) for value in values]


def classify_record(row: list[str], columns: dict[str, int], number: int) -> tuple[int, int]:
    drug_cell = row[columns[DRUG_COLUMN]]
    drugs = parse_ids(drug_cell, DRUG_COLUMN, number)
    require(len(drugs) == 1 and drugs[0] in EXPECTED_COUNTS,
            f"Data record {number}: unexpected {DRUG_COLUMN} cell {drug_cell!r}; "
            "expected exactly one ID in {5, 43, 86}")
    id_cell, fraction_cell = row[columns[SOLVENT_COLUMN]], row[columns[FRACTION_COLUMN]]
    ids = parse_ids(id_cell, SOLVENT_COLUMN, number)
    fractions = parse_list(fraction_cell, FRACTION_COLUMN, number)
    require(len(ids) == len(fractions),
            f"Data record {number}: solvent ID/fraction lengths differ: {id_cell!r}, {fraction_cell!r}")
    require(len(set(ids)) == len(ids),
            f"Data record {number}: repeated solvent IDs in {id_cell!r}")
    require(all(type(value) in (int, float) and 0 <= value <= 1 and math.isfinite(value)
                for value in fractions),
            f"Data record {number}: invalid {FRACTION_COLUMN} cell {fraction_cell!r}; "
            "fractions must be finite numbers in [0, 1]")
    total = math.fsum(fractions)
    require(math.isclose(total, 1.0, rel_tol=0.0, abs_tol=SUM_TOLERANCE),
            f"Data record {number}: fractions sum to {total!r}, outside absolute tolerance "
            f"{SUM_TOLERANCE}: {fraction_cell!r}")
    active = sum(value > 0 for value in fractions)
    require(active in (1, 2),
            f"Data record {number}: unexpected {active} active solvents: {id_cell!r}, {fraction_cell!r}")
    return drugs[0], active


def filter_records(header: list[str], rows: list[list[str]]) -> FilterResult:
    columns = validate_header(header)
    result = FilterResult([], [], [], {cid: {"input": 0, "excluded": 0, "retained": 0}
                                     for cid in EXPECTED_COUNTS})
    for number, row in enumerate(rows, 1):
        require(len(row) == len(header), f"Data record {number}: incorrect field count")
        cid, active = classify_record(row, columns, number)
        result.counts[cid]["input"] += 1
        kind = "retained" if active == 2 else "excluded"
        result.counts[cid][kind] += 1
        if active == 2:
            result.retained.append(row.copy())
            result.retained_positions.append(number)
        else:
            result.excluded_positions.append(number)
    require(len(result.retained) + len(result.excluded_positions) == len(rows),
            "Input/retained/excluded record totals do not balance")
    require(all(count["input"] == count["retained"] + count["excluded"]
                for count in result.counts.values()), "Per-drug record totals do not balance")
    return result


def check_expected_counts(result: FilterResult) -> None:
    actual = {cid: (count["input"], count["excluded"], count["retained"])
              for cid, count in result.counts.items()}
    require(actual == EXPECTED_COUNTS,
            f"Counts disagree with the previous audit (input, excluded, retained): "
            f"actual={actual}; expected={EXPECTED_COUNTS}; output publication refused")


def verify_output(path: Path, header: list[str], expected: FilterResult) -> dict[str, str]:
    actual_header, actual = read_csv(path)
    require(actual_header == header, "Read-back header or column order differs from the input")
    require(actual == expected.retained,
            "Read-back records differ from the complete expected selection in fields, count or order")
    verified = filter_records(actual_header, actual)
    require(not verified.excluded_positions, "Read-back output contains pure-solvent observations")
    require(all(verified.counts[cid]["input"] == expected.counts[cid]["retained"]
                for cid in EXPECTED_COUNTS), "Read-back per-drug counts differ")
    return {"Exactly two active solvents in every output record": "PASS",
            "Only singleton drug IDs 5, 43 and 86": "PASS",
            "Every qualifying input occurrence retained exactly once": "PASS",
            "Complete header and original column order": "PASS",
            "Every retained field and relative record order": "PASS",
            "Read-back total and per-drug counts": "PASS"}


def check_input(path: Path, original_hash: str) -> None:
    require(sha256(path) == original_hash, "Input SHA-256 changed; output publication refused")


def guard_output(output: Path, input_path: Path) -> None:
    require(output.name == OUTPUT_NAME, f"Output filename must be {OUTPUT_NAME}")
    require(not output.is_symlink(), f"Refusing symbolic-link output: {output}")
    require(output.resolve().is_relative_to(REPO_ROOT.resolve()), "Output must stay inside the repository")
    for protected in (input_path, Path(__file__)):
        require(output.resolve() != protected.resolve(), f"Output aliases protected input: {protected}")
        if output.exists() and protected.exists():
            require(not os.path.samefile(output, protected), f"Output aliases protected input: {protected}")
    require(not output.exists() or output.is_file(), f"Output is not a regular file: {output}")


def publish_output(input_path: Path, output: Path, input_hash: str,
                   header: list[str], result: FilterResult) -> tuple[str, dict[str, str]]:
    guard_output(output, input_path)
    check_input(input_path, input_hash)
    if output.exists():
        existing_hash = sha256(output)
        try:
            checks = verify_output(output, header, result)
        except (OSError, ValueError, csv.Error, UnicodeError) as exc:
            raise ValueError(f"Existing output conflicts with the verified result; refusing replacement: {output}: {exc}") from exc
        check_input(input_path, input_hash)
        require(sha256(output) == existing_hash, "Existing output changed during verification")
        return "Existing output exactly matches; reused unchanged", checks

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.parent / f".binary-subset-{uuid.uuid4().hex}.tmp"
    created = published = False
    try:
        with temporary.open("x", encoding="utf-8", newline="") as handle:
            created = True
            writer = csv.writer(handle)
            writer.writerow(header)
            writer.writerows(result.retained)
        verify_output(temporary, header, result)
        check_input(input_path, input_hash)
        guard_output(output, input_path)
        # Same-directory hard-link publication refuses even a racing output.
        os.link(temporary, output)
        published = True
        checks = verify_output(output, header, result)
        check_input(input_path, input_hash)
        return "Published and completely read back new binary-mixture subset", checks
    except BaseException:
        if published:
            output.unlink(missing_ok=True)
        raise
    finally:
        if created:
            temporary.unlink(missing_ok=True)


def self_checks() -> int:
    """Exercise composition rules and field preservation without writing files."""
    header = ["value", DRUG_COLUMN, SOLVENT_COLUMN, FRACTION_COLUMN]
    columns = validate_header(header)
    valid = [("[127.0]", "[1.0]", 1), ("[127, 130]", "[0.0, 1.0]", 1),
             ("[127, 130]", "[1.0, 0.0]", 1), ("[127, 130]", "[0.2, 0.8]", 2),
             ("[127, 130]", "[1e-12, 0.999999999999]", 2),
             ("[127, 130, 145]", "[0.2, 0.8, 0.0]", 2),
             ("[127, 130]", "[0.5, 0.5000005]", 2)]
    checks = 0
    for ids, fractions, active in valid:
        require(classify_record(["-9.000", "[5.0]", ids, fractions], columns, 7) == (5, active),
                "Valid composition self-check failed")
        checks += 1
    invalid = [("", "[1]"), ("127", "[1]"), ("(127,)", "[1]"), ("[]", "[]"),
               ("[[127]]", "[1]"), ("['127']", "[1]"), ("[True]", "[1]"),
               ("[None]", "[1]"), ("[-127]", "[1]"), ("[127.9]", "[1]"),
               ("[1e309]", "[1]"), ("[127,127.0]", "[0.2,0.8]"),
               ("[127,130]", "[1]"), ("[127]", "1"), ("[127]", "(1,)"),
               ("[127]", "[True]"), ("[127]", "['1']"), ("[127]", "[None]"),
               ("[127]", "[[1]]"), ("[127]", "[1e309]"), ("[127]", "[nan]"),
               ("[127]", "[-0.1]"), ("[127]", "[1.1]"), ("[127]", "[0.9]"),
               ("[127,130]", "[0,0]"), ("[127,130]", "[0.5,0.5000015]"),
               ("[127,130,145]", "[0.2,0.3,0.5]"), ("[127", "[1]")]
    for ids, fractions in invalid:
        try:
            classify_record(["0", "[5]", ids, fractions], columns, 7)
        except ValueError as exc:
            require("record 7" in str(exc), "Errors must identify the data record")
        else:
            raise ValueError(f"Invalid composition accepted: {ids!r}, {fractions!r}")
        checks += 1
    for drug in ("", "5", "(5,)", "[]", "[5,43]", "[5,5]", "[15]", "[5.9]",
                 "[True]", "['5']", "[None]", "[-5]", "[1e309]", "[[5]]", "[5"):
        try:
            classify_record(["0", drug, "[127]", "[1]"], columns, 7)
        except ValueError as exc:
            require("record 7" in str(exc), "Drug errors must identify the data record")
        else:
            raise ValueError(f"Invalid drug accepted: {drug!r}")
        checks += 1
    rows = [["-100.000", "[5.0]", "[127,130]", "[0.2,0.8]"],
            ["1.00", "[43]", "[127,130]", "[0,1]"],
            ["2.00", "[86]", "[127,130]", "[0.2,0.8]"],
            ["-100.000", "[5.0]", "[127,130]", "[0.2,0.8]"]]
    snapshot = [row.copy() for row in rows]
    result = filter_records(header, rows)
    require(result.retained == [rows[i] for i in (0, 2, 3)] and rows == snapshot
            and result.retained_positions == [1, 3, 4] and result.excluded_positions == [2],
            "Filtering must preserve every qualifying occurrence, fields and order")
    return checks + 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path, help="Existing three-drug CSV; relative to the repository root")
    parser.add_argument("--output", type=Path, help="Binary-only CSV inside this repository; defaults beside input")
    args = parser.parse_args(argv)
    try:
        print(f"In-memory composition/parser/filter checks: {self_checks()} PASS")
        input_path = locate_input(args.input)
        default_output = (input_path.parent if input_path.is_relative_to(REPO_ROOT.resolve())
                          else REPO_ROOT / "data/processed") / OUTPUT_NAME
        # Retain the unresolved path until guard_output can reject a symlink.
        output = args.output if args.output is not None else default_output
        output = output if output.is_absolute() else REPO_ROOT / output
        guard_output(output, input_path)
        print(f"Resolved input: {input_path}")
        print(f"Resolved output: {output.resolve()}")
        input_hash = sha256(input_path)
        header, rows = read_csv(input_path)
        result = filter_records(header, rows)
        print("Drug ID | Drug | Input | Pure excluded | Binary retained")
        for cid, counts in result.counts.items():
            print(f"{cid} | {DRUG_NAMES[cid]} | {counts['input']} | {counts['excluded']} | {counts['retained']}")
        print(f"Total input: {len(rows)}")
        print(f"Total pure-solvent observations excluded: {len(result.excluded_positions)}")
        print(f"Total binary-mixture observations retained: {len(result.retained)}")
        check_expected_counts(result)
        outcome, checks = publish_output(input_path, output, input_hash, header, result)
        print(outcome)
        print("COMPLETE FINAL-FILE READ-BACK VERIFICATION")
        print("Every excluded input record has exactly one active solvent: PASS")
        print("Overall and per-drug input = retained + excluded: PASS")
        print("Previous audit's expected counts: PASS")
        for label, status in checks.items():
            print(f"{label}: {status}")
        print(f"Original columns preserved: PASS ({len(header)} columns)")
        print(f"Input SHA-256 unchanged: PASS ({input_hash})")
        print(f"Output SHA-256: {sha256(output)}")
        return 0
    except (OSError, ValueError, csv.Error, UnicodeError) as exc:
        print(f"Binary-mixture filtering stopped: {exc}", file=sys.stderr)
        print("Verification: FAIL; no successful output is claimed.", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
