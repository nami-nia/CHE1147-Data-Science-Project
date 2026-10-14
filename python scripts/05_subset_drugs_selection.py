#!/usr/bin/env python3
"""Select all single-drug records for IDs 5, 43 and 86 without cleaning them.

Python 3.10+; standard library only. From the project repository root:

    py -B "python scripts/05_subset_drugs_selection.py"

Defaults are resolved from this script's repository, independent of the current
working directory. Default input discovery requires exactly one file named
id_mismatch_repaired_duplicate_removed_DrugSolubilityData.csv in this repository.
Use --input to choose explicitly if discovery finds multiple copies. All relative
--input, --output, --compounds and --raw-metadata paths are repository-relative.
The default output is data/processed/subset_drugs_DrugSolubilityData.csv.
The default compound/raw metadata belong to the sibling drug-solubility dataset.
With an explicit --compounds, raw metadata is inferred from its dataset directory
unless --raw-metadata is supplied. No metadata is downloaded or modified.

Every drug-ID cell is parsed and validated, including excluded records. Select
only a singleton list containing an exactly integral numeric ID in {5, 43, 86}.
All original field strings, columns and selected-record order are preserved;
no other column is used to filter. Multi-ID lists remain unchanged and excluded.

Lightweight parser and selection checks run in memory before processing. A new
CSV is written to a temporary file in its output directory, read back in full,
and published only after validation and a repeated input SHA-256 check. Existing
outputs are reused only if they exactly match the complete expected selection.
Otherwise replacement requires explicit --overwrite (never implicit).
Output paths must stay inside the repository and cannot alias any source file.
No reports, permanent fixtures, or Git operations are created or performed.
Exit codes: 0 = complete selection verified; 2 = selection or verification failure.
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
from collections import Counter, defaultdict
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
INPUT_NAME = "id_mismatch_repaired_duplicate_removed_DrugSolubilityData.csv"
DEFAULT_OUTPUT = REPO_ROOT / "data/processed/subset_drugs_DrugSolubilityData.csv"
DEFAULT_COMPOUNDS = REPO_ROOT.parent / "chemixhub/datasets/drug-solubility/processed_data/compounds.csv"
EXPECTED_DRUGS = {5: "Itraconazole", 43: "Posaconazole", 86: "Ketoconazole"}
DRUG_COLUMN = "cmp_ids_drug"


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
    if requested is not None:
        path = resolve_path(requested)
        require(path.name == INPUT_NAME, f"Input filename must be {INPUT_NAME}: {path}")
        require(path.is_file(), f"Required cleaned input is missing: {path}")
        return path
    copies = sorted(p.resolve() for p in REPO_ROOT.rglob(INPUT_NAME)
                    if p.is_file() and ".git" not in p.relative_to(REPO_ROOT).parts)
    require(bool(copies), f"Required cleaned input {INPUT_NAME} is missing in {REPO_ROOT}")
    require(len(copies) == 1, "Multiple cleaned input copies found; choose with --input:\n"
            + "\n".join(str(p) for p in copies))
    return copies[0]


def read_csv(path: Path) -> tuple[list[str], list[list[str]]]:
    """Keep decoded cells as strings; reject ambiguous headers and bad widths."""
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, strict=True)
        header = next(reader, [])
        require(bool(header) and all(name.strip() for name in header),
                f"{path}: empty or unnamed CSV columns")
        require(len(header) == len(set(header)), f"{path}: duplicate column names")
        rows = []
        try:
            for record_number, row in enumerate(reader, 1):
                require(len(row) == len(header),
                        f"{path}: data record {record_number} has {len(row)} fields; expected {len(header)}")
                rows.append(row)
        except csv.Error as exc:
            raise ValueError(f"{path}: invalid CSV quoting at data record {len(rows) + 1}: {exc}") from exc
    return header, rows


def numeric_id(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("IDs must be numeric integers or integral finite floats")
    if isinstance(value, float) and (not math.isfinite(value) or not value.is_integer()):
        raise ValueError("IDs must be finite and exactly integer-valued")
    if value < 0:
        raise ValueError("IDs must be non-negative")
    return int(value)


def parse_drug_ids(text: str, record_number: int) -> list[int]:
    try:
        values = ast.literal_eval(text)
        if not isinstance(values, list) or not values:
            raise ValueError("expected a nonempty list")
        return [numeric_id(value) for value in values]
    except (ValueError, SyntaxError, TypeError, RecursionError, MemoryError, OverflowError) as exc:
        raise ValueError(f"Input data record {record_number}: invalid {DRUG_COLUMN} cell {text!r}: {exc}") from exc


def select_records(header: list[str], rows: list[list[str]]) -> tuple[list[list[str]], dict[int, int]]:
    require(len(header) == len(set(header)) and header.count(DRUG_COLUMN) == 1,
            f"Input header must contain {DRUG_COLUMN} exactly once and have unique columns")
    column = header.index(DRUG_COLUMN)
    selected = []
    counts = {cid: 0 for cid in EXPECTED_DRUGS}
    for number, row in enumerate(rows, 1):
        require(len(row) == len(header), f"Input data record {number}: incorrect field count")
        ids = parse_drug_ids(row[column], number)
        if len(ids) == 1 and ids[0] in EXPECTED_DRUGS:
            selected.append(row.copy())
            counts[ids[0]] += 1
    return selected, counts


def verify_mappings(compound_header: list[str], compounds: list[list[str]],
                    raw_header: list[str], raw_rows: list[list[str]]) -> dict[int, tuple[str, str]]:
    """Verify explicit IDs and names; exact SMILES provides a fallback/cross-check."""
    by_id = defaultdict(list)
    if "compound_id" in compound_header:
        column = compound_header.index("compound_id")
        for number, row in enumerate(compounds, 1):
            try:
                cid = numeric_id(ast.literal_eval(row[column].strip()))
            except (ValueError, SyntaxError, TypeError, OverflowError) as exc:
                raise ValueError(f"Compound metadata record {number}: invalid explicit compound_id {row[column]!r}") from exc
            by_id[cid].append(dict(zip(compound_header, row)))
    raw_names = defaultdict(set)
    if "Drug" in raw_header and "SMILES_0" in raw_header:
        name_column, smiles_column = raw_header.index("Drug"), raw_header.index("SMILES_0")
        for row in raw_rows:
            if row[name_column].strip() and row[smiles_column].strip():
                raw_names[row[smiles_column]].add(row[name_column].strip().casefold())
    results = {}
    for cid, expected in EXPECTED_DRUGS.items():
        entries = by_id.get(cid, [])
        if not entries:
            results[cid] = ("UNVERIFIED", "no unambiguous explicit-ID metadata available")
            continue
        names = {r.get("name", "").strip().casefold() for r in entries if r.get("name", "").strip()}
        smiles = {r.get("smiles", "") for r in entries if r.get("smiles", "").strip()}
        if len(names) > 1 or len(smiles) > 1 or (names and names != {expected.casefold()}):
            results[cid] = ("CONFLICT", f"metadata disagrees with {expected}: names={sorted(names)}; distinct SMILES={len(smiles)}")
            continue
        matched_names = raw_names.get(next(iter(smiles)), set()) if smiles else set()
        if matched_names and matched_names != {expected.casefold()}:
            results[cid] = ("CONFLICT", f"exact SMILES matches conflicting raw drug names: {sorted(matched_names)}")
        elif names:
            results[cid] = ("VERIFIED", "explicit compound_id/name" + ("; exact raw SMILES/name cross-check" if matched_names else ""))
        elif matched_names:
            results[cid] = ("VERIFIED", "explicit compound_id with exact, unambiguous raw SMILES/name match")
        else:
            results[cid] = ("UNVERIFIED", "available ID/SMILES metadata cannot independently establish the drug name")
    return results


def load_metadata(compounds_path: Path, raw_path: Path) -> dict[int, tuple[str, str]]:
    compound_header, compounds = read_csv(compounds_path) if compounds_path.is_file() else ([], [])
    raw_header, raw_rows = read_csv(raw_path) if raw_path.is_file() else ([], [])
    return verify_mappings(compound_header, compounds, raw_header, raw_rows)


def guard_output(output: Path, protected: list[Path]) -> None:
    resolved = output.resolve()
    require(resolved.is_relative_to(REPO_ROOT.resolve()), f"Output must remain inside the project repository: {resolved}")
    for path in protected:
        require(resolved != path.resolve(), f"Output would overwrite a protected source: {path}")
        if output.exists() and path.exists():
            require(not os.path.samefile(output, path), f"Output aliases a protected source: {path}")
    require(not output.is_symlink(), f"Refusing a symbolic-link output: {output}")
    require(not output.exists() or output.is_file(), f"Output is not a regular file: {output}")


def verify_output(path: Path, header: list[str], expected: list[list[str]],
                  counts: dict[int, int]) -> dict[str, str]:
    actual_header, actual = read_csv(path)
    require(actual_header == header, "Read-back output header or column order changed")
    require(actual == expected, "Read-back output does not match the complete selected input, field by field and in order")
    column = header.index(DRUG_COLUMN)
    observed = Counter()
    for number, row in enumerate(actual, 1):
        ids = parse_drug_ids(row[column], number)
        require(len(ids) == 1 and ids[0] in EXPECTED_DRUGS,
                f"Output data record {number} violates the single-target-ID rule")
        observed[ids[0]] += 1
    require(len(actual) == len(expected) == sum(counts.values()), "Selected counts do not add up")
    require(all(counts[cid] > 0 and observed[cid] == counts[cid] for cid in EXPECTED_DRUGS),
            "Drug coverage or per-ID counts failed")
    return {"Selection completeness": "PASS", "Selection exclusivity": "PASS",
            "Total and per-ID counts": "PASS", "All three requested drugs present": "PASS",
            "Original header and column order": "PASS", "All original field contents": "PASS",
            "Relative order of selected records": "PASS"}


def check_input(path: Path, original_hash: str) -> None:
    require(sha256(path) == original_hash, "Input SHA-256 changed during selection; output publication refused")


def publish_subset(input_path: Path, output: Path, protected: list[Path], input_hash: str,
                   header: list[str], selected: list[list[str]], counts: dict[int, int],
                   overwrite: bool) -> tuple[str, dict[str, str]]:
    guard_output(output, protected)
    previous_hash = sha256(output) if output.is_file() else None
    if previous_hash is not None:
        try:
            checks = verify_output(output, header, selected, counts)
        except (OSError, ValueError, csv.Error, UnicodeError):
            require(overwrite, "Existing subset fails complete verification; refusing replacement without explicit --overwrite")
        else:
            check_input(input_path, input_hash)
            require(sha256(output) == previous_hash, "Existing subset changed during verification")
            return "Reused existing, fully verified subset; no file rewritten", checks
    check_input(input_path, input_hash)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.parent / f".subset-selection-{uuid.uuid4().hex}.tmp"
    backup = output.parent / f".subset-backup-{uuid.uuid4().hex}.tmp"
    temporary_created = backup_created = published = False
    try:
        with temporary.open("x", encoding="utf-8", newline="") as handle:
            temporary_created = True
            writer = csv.writer(handle)
            writer.writerow(header)
            writer.writerows(selected)
        verify_output(temporary, header, selected, counts)
        check_input(input_path, input_hash)
        guard_output(output, protected)
        if previous_hash is not None:
            require(overwrite and output.is_file() and sha256(output) == previous_hash,
                    "Existing output changed before publication; refusing replacement")
            with backup.open("xb") as handle:
                backup_created = True
                handle.write(output.read_bytes())
            require(sha256(backup) == previous_hash and sha256(output) == previous_hash,
                    "Existing output changed while preparing safe replacement")
            os.replace(temporary, output)
        else:
            # Atomic publication without overwriting a file created in a race.
            # Both paths are in the same output directory/filesystem.
            os.link(temporary, output)
        published = True
        checks = verify_output(output, header, selected, counts)
        check_input(input_path, input_hash)
        return "Published new subset" if previous_hash is None else "Replaced subset with explicit --overwrite", checks
    except BaseException:
        if published:
            if backup_created:
                os.replace(backup, output)
            else:
                output.unlink(missing_ok=True)
        raise
    finally:
        if temporary_created:
            temporary.unlink(missing_ok=True)
        if backup_created:
            backup.unlink(missing_ok=True)


def self_checks() -> int:
    """Small in-memory parser/selection checks; no files, libraries or fixtures."""
    checks = 0
    for text, expected in [("[5]", [5]), ("[43]", [43]), ("[86]", [86]),
                           ("[5.0]", [5]), ("[15]", [15]), ("[5, 43]", [5, 43]),
                           ("[5, 5]", [5, 5]), ("[0]", [0])]:
        require(parse_drug_ids(text, 7) == expected, f"Parser self-check failed: {text}")
        checks += 1
    invalid = ["", " ", "not a list", "5", "(5,)", "[]", "[[5]]", "['5']", "[True]",
               "[False]", "[None]", "[-5]", "[5.9]", "[1e309]", "[nan]", "[5", "[5, 5.9]", "[5+0j]"]
    for text in invalid:
        try:
            parse_drug_ids(text, 7)
        except ValueError as exc:
            require("record 7" in str(exc) and repr(text) in str(exc), "Parser errors must identify the record and offending value")
        else:
            raise ValueError(f"Invalid cell unexpectedly accepted: {text!r}")
        checks += 1
    header = ["value", DRUG_COLUMN, "note"]
    rows = [["-100.000", "[5.0]", 'comma, quote " and newline\nretained'],
            ["99", "[15]", ""], ["0", "[5, 43]", ""], ["0", "[5, 5]", ""],
            ["1.00", "[43]", "pure endpoint"], ["2", "[86]", ""],
            ["-100.000", "[5.0]", 'comma, quote " and newline\nretained']]
    snapshot = [r.copy() for r in rows]
    selected, counts = select_records(header, rows)
    require(selected == [rows[i] for i in (0, 4, 5, 6)] and counts == {5: 2, 43: 1, 86: 1}
            and rows == snapshot, "Selection must preserve duplicates, extreme values, original strings and order")
    checks += 1
    try:
        select_records(header, [["1", "[15]", ""], ["2", "[True]", ""]])
    except ValueError as exc:
        require("record 2" in str(exc), "Excluded records must also be parsed")
    else:
        raise ValueError("Malformed excluded cell was silently skipped")
    return checks + 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path, help="Exact cleaned CSV; relative paths start at the repository root")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Subset CSV inside this repository")
    parser.add_argument("--compounds", type=Path, default=DEFAULT_COMPOUNDS, help="Associated compound-ID/name/SMILES metadata")
    parser.add_argument("--raw-metadata", type=Path, help="Associated raw Drug/SMILES_0 metadata for exact-name verification")
    parser.add_argument("--overwrite", action="store_true", help="Explicitly replace an existing subset only after all checks")
    args = parser.parse_args(argv)
    try:
        print(f"In-memory parser and selection checks: {self_checks()} PASS")
        input_path = locate_input(args.input)
        output = resolve_path(args.output)
        compounds = resolve_path(args.compounds)
        raw = resolve_path(args.raw_metadata) if args.raw_metadata else compounds.parent.parent / "raw_data/DrugSolubilityData.csv"
        protected = [input_path, compounds, raw, Path(__file__).resolve()]
        guard_output(output, protected)
        print(f"Resolved input: {input_path}")
        print(f"Resolved output: {output}")
        print(f"Compound metadata: {compounds}")
        print(f"Raw drug metadata: {raw}")
        mappings = load_metadata(compounds, raw)
        for cid, (status, detail) in mappings.items():
            print(f"ID {cid} -> {EXPECTED_DRUGS[cid]}: {status}; {detail}")
            if status == "UNVERIFIED":
                print(f"WARNING: ID {cid}/{EXPECTED_DRUGS[cid]} name mapping could not be independently verified; selection uses the requested ID.")
        require(not any(status == "CONFLICT" for status, detail in mappings.values()),
                "Contradictory/conflicting drug metadata; no subset will be created or replaced")
        input_hash = sha256(input_path)
        header, rows = read_csv(input_path)
        selected, counts = select_records(header, rows)
        print(f"Input records: {len(rows):,}")
        for cid in EXPECTED_DRUGS:
            print(f"Selected ID {cid} ({EXPECTED_DRUGS[cid]}): {counts[cid]:,}")
        print(f"Total selected: {len(selected):,}")
        print(f"Total excluded: {len(rows) - len(selected):,}")
        print(f"Unique subset IDs: {sorted(cid for cid, count in counts.items() if count)}")
        require(all(counts[cid] > 0 for cid in EXPECTED_DRUGS),
                "Requested three-drug selection is incomplete: at least one requested drug has zero observations; no final CSV created or replaced")
        outcome, checks = publish_subset(input_path, output, protected, input_hash, header,
                                         selected, counts, args.overwrite)
        print(outcome)
        print("COMPLETE FINAL-FILE READ-BACK VERIFICATION")
        for label, status in checks.items():
            print(f"{label}: {status}")
        print("Input SHA-256 unchanged: PASS")
        print(f"Input SHA-256: {input_hash}")
        print(f"Output SHA-256: {sha256(output)}")
        return 0
    except (OSError, ValueError, csv.Error, UnicodeError) as exc:
        print(f"Drug selection stopped: {exc}", file=sys.stderr)
        print("Verification: FAIL; no successful selection is claimed.", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
