#!/usr/bin/env python3
"""Read-only Step 2 validation for the CHE1147 drug-solubility dataset.

Python 3.10+; standard library only. No installation of extra packages is needed.
Run from the GROUP repository root, not from the upstream dataset repository:

    python scripts/02_validate_values.py --dataset-dir chemixhub/datasets/drug-solubility

Add --report-dir reports/data_cleaning/02_validation to write validation reports.
The three CSV inputs and croissant.json are never modified. No records are
repaired, removed, averaged, reordered, or converted to another concentration unit.

Exit codes: 0 = no structural errors (human-review items may remain);
            1 = validation completed and found structural errors;
            2 = input, schema, or output-path error.

Scope: preserve the supplied column meanings and test separate drug/solvent
normalization. The README requires corresponding ID/fraction order. The supplied
processing script removes zero fractions and builds lists independently; these
checks detect resulting mismatches. Logarithm and unit semantics still require
source review. This script does not validate chemical structures with RDKit.
"""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

PAIR_TOL = 1e-12  # Only accommodates floating-point serialization differences.
ISSUE_FIELDS = ["severity", "file", "data_row", "raw_label", "check", "detail"]
DUPLICATE_FIELDS = [
    "group", "raw_data_row", "raw_label", "drug", "temperature_K", "LogS",
    "solvent_id_fraction_pairs",
]


def load_csv(path: Path, required: set[str]) -> list[dict]:
    """Read a CSV and reject an unexpected schema or malformed record width."""
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        header = reader.fieldnames or []
        if len(header) != len(set(header)):
            raise ValueError(f"{path.name}: duplicate column names")
        missing = required - set(header)
        if missing:
            raise ValueError(f"{path.name}: missing columns {sorted(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{path.name}: no data records")
    if any(None in row or None in row.values() for row in rows):
        raise ValueError(f"{path.name}: a record has the wrong number of fields")
    return rows


def number(value: object) -> float | None:
    """Return a finite number, or None for invalid/missing/nonfinite input."""
    if isinstance(value, bool):
        return None
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError, OverflowError):
        return None


def canonical_pairs(ids: list, fractions: list) -> tuple:
    """Compare ID/fraction pairs together; ignore only exactly zero fractions.

    Call only after equal-length and numeric-validity checks. Do not normalize
    fractions, round nearby compositions, or infer missing fractions.
    """
    return tuple(sorted((int(i), float(x)) for i, x in zip(ids, fractions) if x != 0))


def same_pairs(left: tuple, right: tuple) -> bool:
    return len(left) == len(right) and all(
        i == j and math.isclose(x, y, rel_tol=0, abs_tol=PAIR_TOL)
        for (i, x), (j, y) in zip(left, right)
    )


def validate(raw_path: Path, processed_path: Path, compounds_path: Path,
             metadata_path: Path, sum_tol: float) -> tuple[str, list, list]:
    """Run checks in memory. Return summary text, issues, duplicate candidates."""
    raw = load_csv(raw_path, {
        "LogS", "Temperature, K", "Drug", "SMILES_0", "SMILES_1", "SMILES_2",
        "comp_0", "comp_1", "comp_2",
    })
    processed = load_csv(processed_path, {
        "value", "Temperature, K", "property", "unit", "cmp_ids_solvent",
        "cmp_mole_fractions_solvent", "cmp_ids_drug", "cmp_mole_fractions_drug",
    })
    compounds = load_csv(compounds_path, {"compound_id", "smiles"})
    with metadata_path.open(encoding="utf-8-sig") as handle:
        metadata = json.load(handle)
    if not isinstance(metadata, dict):
        raise ValueError("croissant.json: expected a JSON object")

    # Do not silently choose between ambiguous ID/SMILES mappings.
    smiles_to_id = {}
    known_ids = set()
    for row in compounds:
        cid = number(row["compound_id"])
        smi = row["smiles"].strip()
        if cid is None or cid < 0 or not cid.is_integer() or not smi:
            raise ValueError("compounds.csv: invalid ID or empty SMILES")
        cid = int(cid)
        if cid in known_ids or smi in smiles_to_id:
            raise ValueError("compounds.csv: duplicate ID or ambiguous SMILES mapping")
        known_ids.add(cid)
        smiles_to_id[smi] = cid

    issues = []

    def add(check, detail, row="", label="", file="processed", severity="ERROR"):
        issues.append(dict(zip(ISSUE_FIELDS, [severity, file, row, label, check, detail])))

    def parse_list(text, field, row, label, is_ids=False):
        try:
            if len(text) > 10000:
                raise ValueError("list cell is unexpectedly large")
            values = ast.literal_eval(text)
        except (ValueError, SyntaxError, TypeError, RecursionError):
            add("INVALID_LIST", f"{field}: cannot parse {text[:160]!r}", row, label)
            return None
        if not isinstance(values, list) or not values:
            add("INVALID_LIST", f"{field}: expected a nonempty list", row, label)
            return None
        if any(isinstance(x, bool) or not isinstance(x, (int, float))
               or number(x) is None for x in values):
            add("INVALID_LIST_ITEM", f"{field}: invalid numeric item", row, label)
            return None
        if is_ids:
            if any(x < 0 or x != int(x) or int(x) not in known_ids for x in values):
                add("INVALID_ID", f"{field}: noninteger, negative, or unknown ID", row, label)
                return None
            if len(set(values)) != len(values):
                add("REPEATED_ID", f"{field}: repeated ID within list", row, label)
                return None
        elif (any(x < 0 or x > 1 for x in values)
              or not math.isclose(sum(values), 1, rel_tol=0, abs_tol=sum_tol)):
            add("INVALID_FRACTIONS", f"{field}: out of bounds or sum={sum(values):.12g}",
                row, label)
            return None
        return values

    # Raw compositions provide an independent reference for the processed lists.
    expected = []
    condition_groups = defaultdict(list)
    for row_n, row in enumerate(raw, start=1):
        label = row.get("Label", "")
        comps = [number(row[f"comp_{k}"]) for k in range(3)]
        y, temperature = number(row["LogS"]), number(row["Temperature, K"])
        if y is None or temperature is None or temperature <= 0:
            add("INVALID_RAW_NUMERIC", "Invalid LogS or temperature", row_n, label, "raw")
        valid = all(x is not None and 0 <= x <= 100 for x in comps)
        if valid:
            valid = (math.isclose(comps[0], 100, abs_tol=100 * sum_tol, rel_tol=0)
                     and math.isclose(comps[1] + comps[2], 100,
                                      abs_tol=100 * sum_tol, rel_tol=0))
        if not valid:
            add("INVALID_RAW_COMPOSITION", "Expected separate drug/solvent totals of 100",
                row_n, label, "raw")
            expected.append(None)
            continue
        pairs = []
        for k, percent in enumerate(comps):
            if percent == 0:
                pairs.append(None)  # A missing SMILES is allowed for an absent component.
                continue
            cid = smiles_to_id.get(row[f"SMILES_{k}"].strip())
            if cid is None:
                add("RAW_COMPOUND_NOT_MAPPED", f"Positive comp_{k} has no SMILES mapping",
                    row_n, label, "raw")
                valid = False
            pairs.append((cid, percent / 100))
        if not valid:
            expected.append(None)
            continue
        drug_pairs = (pairs[0],)
        solvent_pairs = tuple(sorted(pair for pair in pairs[1:] if pair is not None))
        expected.append({"drug": drug_pairs, "solvent": solvent_pairs})
        if y is not None and temperature is not None and temperature > 0:
            # Exact condition matching, not an outlier test or duplicate deletion rule.
            condition_groups[(drug_pairs, temperature, solvent_pairs)].append(row_n)

    aligned_count = len(raw) == len(processed)
    if not aligned_count:
        add("ROW_COUNT_MISMATCH", "Raw/processed counts differ; positional comparison skipped")
    copied_targets = copied_temperatures = 0
    for row_n, row in enumerate(processed, start=1):
        source = raw[row_n - 1] if aligned_count else None
        label = source.get("Label", "") if source is not None else ""
        y = number(row["value"])
        temperature = number(row["Temperature, K"])
        if y is None:
            add("INVALID_VALUE", "value is missing, nonnumeric, or nonfinite", row_n, label)
        if temperature is None or temperature <= 0:
            add("INVALID_TEMPERATURE", "Temperature must be finite and above 0 K", row_n, label)
        if row["property"].strip() != "Log solubility":
            add("UNEXPECTED_PROPERTY", repr(row["property"]), row_n, label)
        if not row["unit"].strip():
            add("EMPTY_UNIT", "Unit label is empty", row_n, label)
        if source is not None:
            for key, value, raw_key in [("TARGET", y, "LogS"),
                                        ("TEMPERATURE", temperature, "Temperature, K")]:
                raw_value = number(source[raw_key])
                if value is not None and raw_value is not None and value == raw_value:
                    if key == "TARGET":
                        copied_targets += 1
                    else:
                        copied_temperatures += 1
                else:
                    add(f"{key}_NOT_PRESERVED", "Check file versions, row order, and values",
                        row_n, label)
        for kind in ("solvent", "drug"):
            id_col = f"cmp_ids_{kind}"
            frac_col = f"cmp_mole_fractions_{kind}"
            ids = parse_list(row[id_col], id_col, row_n, label, is_ids=True)
            fractions = parse_list(row[frac_col], frac_col, row_n, label)
            if ids is None or fractions is None:
                continue
            if len(ids) != len(fractions):
                add(f"{kind.upper()}_LENGTH_MISMATCH",
                    f"{len(ids)} IDs vs {len(fractions)} fractions; "
                    f"IDs={ids}; fractions={fractions}", row_n, label)
                continue  # Never guess the association by zipping unequal lists.
            if source is not None and expected[row_n - 1] is not None:
                actual = canonical_pairs(ids, fractions)
                if not same_pairs(actual, expected[row_n - 1][kind]):
                    add(f"{kind.upper()}_PAIR_MISMATCH",
                        "Processed ID/fraction pairs differ from the raw record; "
                        "check file versions and row order", row_n, label)

    units = sorted({row["unit"].strip() for row in processed})
    description = metadata.get("description", "")
    add("UNIT_SOURCE_REVIEW", f"Observed unit labels: {units}. "
        "Logarithm/concentration definitions require source review; no conversion performed. "
        f"Croissant description: {description}", severity="REVIEW", file="metadata")

    candidates = []
    duplicate_group_count = 0
    for (_, temperature, solvent_pairs), members in condition_groups.items():
        if len(members) < 2:
            continue
        duplicate_group_count += 1
        for row_n in members:
            row = raw[row_n - 1]
            candidates.append(dict(zip(DUPLICATE_FIELDS, [
                duplicate_group_count, row_n, row.get("Label", ""), row["Drug"],
                temperature, row["LogS"], json.dumps(solvent_pairs),
            ])))
    if candidates:
        add("REPEATED_CONDITIONS_REVIEW",
            f"{duplicate_group_count} groups / {len(candidates)} raw rows share drug, "
            "temperature, and nonzero solvent composition. Review, do not auto-delete.",
            severity="REVIEW", file="raw")

    errors = [issue for issue in issues if issue["severity"] == "ERROR"]
    counts = Counter(issue["check"] for issue in errors)
    lines = [
        "STEP 2 - DRUG SOLUBILITY VALIDATION",
        "Mode: read-only (no source data changed)",
        f"Raw records: {len(raw):,}",
        f"Processed records: {len(processed):,}",
        f"Compound records: {len(compounds):,}",
        f"Unchanged target values: {copied_targets:,} / {len(processed):,}",
        f"Unchanged temperatures: {copied_temperatures:,} / {len(processed):,}",
        f"Solvent ID/fraction length mismatches: {counts['SOLVENT_LENGTH_MISMATCH']:,}",
        f"Repeated-condition candidate groups: {duplicate_group_count:,}",
        f"Rows in those groups: {len(candidates):,}",
        f"Structural error findings: {len(errors):,}",
        f"Human-review items: {len(issues) - len(errors):,}",
        "Status: ACTION REQUIRED" if errors else "Status: STRUCTURAL CHECKS PASS; SOURCE REVIEW STILL REQUIRED",
        "", "ERROR COUNTS",
    ]
    lines.extend(f"{key}: {count:,}" for key, count in sorted(counts.items()))
    if not errors:
        lines.append("None")
    lines += ["", "HUMAN REVIEW (not resolved automatically)"]
    lines.extend(f"{i['check']}: {i['detail']}" for i in issues if i["severity"] == "REVIEW")
    lines += [
        "", "METHOD",
        "data_row = 1 means the first CSV record after the header; raw_label is separate.",
        "Raw/processed comparisons use input order, with target, temperature and pair checks.",
        "A comparison failure does not distinguish reordered rows from corrupted values.",
        f"Fraction sum absolute tolerance: {sum_tol:g}; relative tolerance: 0.",
        f"Raw/processed pair tolerance: {PAIR_TOL:g}; target/temperature comparison: exact.",
        "Candidate groups use exact numeric conditions from raw records, ignoring zero components.",
        "Known IDs are checked; chemical structure validity and experimental truth are not.",
        "No logarithm, renormalization, unit correction, repair, or deduplication is performed.",
        "", f"Python: {sys.version.split()[0]}",
        f"Run UTC: {datetime.now(timezone.utc).isoformat(timespec='seconds')}",
        "INPUT SHA256 (identifies the exact input copies)",
    ]
    for path in (raw_path, processed_path, compounds_path, metadata_path):
        lines.append(f"{path.name}: {hashlib.sha256(path.read_bytes()).hexdigest()}")
    return "\n".join(lines) + "\n", issues, candidates


def write_reports(folder: Path, summary: str, issues: list, candidates: list,
                  input_paths: list[Path]) -> None:
    outputs = [folder / name for name in [
        "02_validation_summary.txt", "02_validation_issues.csv", "02_duplicate_candidates.csv",
    ]]
    if any(path.resolve() in {p.resolve() for p in input_paths} for path in outputs):
        raise ValueError("A report path would overwrite an input file")
    folder.mkdir(parents=True, exist_ok=True)
    outputs[0].write_text(summary, encoding="utf-8")
    for path, fields, rows in [(outputs[1], ISSUE_FIELDS, issues),
                               (outputs[2], DUPLICATE_FIELDS, candidates)]:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dataset-dir", type=Path,
                        default=Path("chemixhub/datasets/drug-solubility"),
                        help="Folder containing raw_data/ and processed_data/; relative to the terminal")
    parser.add_argument("--report-dir", type=Path,
                        help="Optional output folder. Omit for console-only validation.")
    parser.add_argument("--sum-tolerance", type=float, default=1e-6)
    args = parser.parse_args()
    if not math.isfinite(args.sum_tolerance) or not 0 < args.sum_tolerance < 1:
        parser.error("--sum-tolerance must be finite and between 0 and 1")
    root = args.dataset_dir
    paths = [root / "raw_data" / "DrugSolubilityData.csv",
             root / "processed_data" / "processed_DrugSolubilityData.csv",
             root / "processed_data" / "compounds.csv",
             root / "processed_data" / "croissant.json"]
    try:
        summary, issues, candidates = validate(*paths, args.sum_tolerance)
        print(summary)
        examples = [i for i in issues if i["severity"] == "ERROR"][:5]
        if examples:
            print("FIRST FIVE STRUCTURAL FINDINGS")
            for issue in examples:
                print(f"{issue['file']} data row {issue['data_row']} "
                      f"(raw Label {issue['raw_label']}): {issue['check']} - {issue['detail']}")
        if args.report_dir:
            write_reports(args.report_dir, summary, issues, candidates, paths)
            print(f"\nReports written to: {args.report_dir.resolve()}")
        else:
            print("\nConsole-only run: no reports written.")
        return 1 if any(i["severity"] == "ERROR" for i in issues) else 0
    except (OSError, ValueError, csv.Error) as exc:
        print(f"Input/report error: {exc}", file=sys.stderr)
        print("Check --dataset-dir and the required CSV column names.", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
