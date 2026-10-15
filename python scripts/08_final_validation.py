#!/usr/bin/env python3
"""Validate the existing analysis population and document the finished review.

Python 3.10+; standard library only. Run from the repository root:

    py -B "python scripts/08_final_validation.py"

Relative --input, --report-dir and --source-dir paths start at this script's
repository, independently of the working directory. The original three-drug
subset is discovered from INPUT_NAME in 06_remove_pure_solvents.py using AST;
it is never guessed. Existing modules are inspected, never executed. This
validator does not recompute quartiles, run IQR screening, filter or write data.
It reconciles existing flags, counts, stored fences and the authoritative review
decision. Duplicate keys exclude all added audit columns; the second comparison
sorts analysis-only solvent ID/fraction pairs with exact numeric equality.

Only reports/validation/final_validation_checks.csv and final_cleaning_summary.md
are written. Identical existing reports are reused; conflicting reports are
refused. Pre-existing datasets, scripts, reports, figures and source files are
hashed before/after. No Git mutations or network requests occur.

PASS/FAIL refer to validation checks. UNRESOLVED identifies documentation or
experimental-origin questions in the historical check record, not a passed
measurement-definition check. The approved source review (2026-10-14) is recorded
separately in reporting text below; reruns preserve its dates, versions and
qualifications without fetching external sources or reading an output template.
SKIPPED identifies a check blocked by a failed prerequisite. Exit 0 means no
failed integrity checks, with any UNRESOLVED questions explicitly retained;
exit 2 means final validation has NOT passed or outputs could not be verified.
"""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import io
import json
import math
import os
import subprocess
import sys
import uuid
from collections import Counter, defaultdict
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation, localcontext
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
REPOSITORY = "nami-nia/CHE1147-Data-Science-Project"
BRANCH = "chenyue-work"
DEFAULT_SOURCE = REPO_ROOT.parent / "chemixhub/datasets/drug-solubility"
DEFAULT_INPUT = REPO_ROOT / "data/processed/subset_drugs_binary_only_DrugSolubilityData.csv"
COLUMNS = ["value", "Temperature, K", "property", "unit", "cmp_ids_solvent",
           "cmp_mole_fractions_solvent", "cmp_ids_drug", "cmp_mole_fractions_drug"]
EXPECTED_COUNTS = {5: 360, 43: 162, 86: 90}
EXPECTED_ORIGINAL = {5: 390, 43: 189, 86: 105}
EXPECTED_EXCLUDED = {5: 30, 43: 27, 86: 15}
DRUG_NAMES = {5: "Itraconazole", 43: "Posaconazole", 86: "Ketoconazole"}
FLAG_COLUMNS = [f"{variable}_iqr_{direction}_by_{level}"
                for level in ("drug", "system") for variable in ("value", "temperature")
                for direction in ("low", "high", "flag")]
CHECK_COLUMNS = ["check_name", "expected_result", "observed_result", "status", "notes"]
COMMAND = 'py -B "python scripts/08_final_validation.py"'


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def show(value: object) -> str:
    return value if isinstance(value, str) else json.dumps(value, sort_keys=True, ensure_ascii=False)


def relative(path: Path) -> str:
    return Path(os.path.relpath(path, REPO_ROOT)).as_posix()


def resolve(path: Path) -> Path:
    return (path if path.is_absolute() else REPO_ROOT / path).resolve()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git(*arguments: str) -> str:
    return subprocess.check_output(["git", "-C", str(REPO_ROOT), *arguments], text=True).strip()


def confirm_repository() -> dict[str, str]:
    root, branch = Path(git("rev-parse", "--show-toplevel")).resolve(), git("branch", "--show-current")
    fetch = git("remote", "get-url", "--all", "origin")
    destination = fetch.removesuffix(".git").rstrip("/")
    permitted = {f"https://github.com/{REPOSITORY}", f"git@github.com:{REPOSITORY}"}
    require(root == REPO_ROOT.resolve() and branch == BRANCH and destination in permitted,
            "Repository root, GitHub origin or branch differs from the approved destination")
    return {"repository": REPOSITORY, "branch": branch, "head": git("rev-parse", "HEAD"),
            "staged": git("diff", "--cached", "--name-status")}


class Checks:
    def __init__(self) -> None:
        self.rows: list[dict[str, str]] = []

    def add(self, name: str, expected: object, observed: object, status: str, notes: str = "") -> None:
        self.rows.append(dict(zip(CHECK_COLUMNS, (name, show(expected), show(observed), status, notes))))

    def compare(self, name: str, expected: object, observed: object, notes: str = "") -> bool:
        passed = observed == expected
        self.add(name, expected, observed, "PASS" if passed else "FAIL", notes)
        return passed

    def attempt(self, name: str, expected: object, function, notes: str = ""):
        try:
            result = function()
        except (OSError, ValueError, KeyError, TypeError, InvalidOperation, SyntaxError) as exc:
            detail = str(exc).replace(str(REPO_ROOT), ".").replace(str(REPO_ROOT.parent), "..")
            self.add(name, expected, detail, "FAIL", notes)
            return None
        self.add(name, expected, "verified", "PASS", notes)
        return result

    def unresolved(self, name: str, observed: str, notes: str) -> None:
        self.add(name, "established by supporting source documentation", observed, "UNRESOLVED", notes)


@dataclass
class Table:
    path: Path
    header: list[str]
    rows: list[list[str]]

    def dictionaries(self) -> list[dict[str, str]]:
        return [dict(zip(self.header, row)) for row in self.rows]


@dataclass
class Observation:
    index: int
    row: list[str]
    drug: int
    pairs: tuple[tuple[int, Decimal], ...]
    value: Decimal
    temperature: Decimal

    @property
    def active(self) -> int:
        return sum(fraction > 0 for _, fraction in self.pairs)

    @property
    def system(self) -> str:
        return "|".join(str(cid) for cid, _ in sorted(self.pairs))


def read_csv(path: Path) -> Table:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, strict=True)
        header = next(reader, [])
        require(bool(header) and all(name.strip() for name in header), f"{relative(path)}: empty/unnamed header")
        require(len(header) == len(set(header)), f"{relative(path)}: ambiguous duplicate column names")
        rows = []
        try:
            for number, row in enumerate(reader, 1):
                require(len(row) == len(header), f"{relative(path)}: data record {number} has incorrect field count")
                rows.append(row)
        except csv.Error as exc:
            raise ValueError(f"{relative(path)}: invalid CSV quoting at data record {len(rows) + 1}: {exc}") from exc
    return Table(path, header, rows)


def parse_list(cell: str, field: str, number: int) -> list:
    try:
        values = ast.literal_eval(cell)
        require(type(values) is list and bool(values), "expected a nonempty list")
        return values
    except (ValueError, SyntaxError, TypeError, OverflowError, RecursionError, MemoryError) as exc:
        raise ValueError(f"Data record {number}: invalid {field}: {cell!r}: {exc}") from exc


def ids(cell: str, field: str, number: int) -> list[int]:
    values = parse_list(cell, field, number)
    require(all(type(v) in (int, float) and v >= 0 and
                (type(v) is int or (math.isfinite(v) and v.is_integer())) for v in values),
            f"Data record {number}: invalid numeric IDs in {field}: {cell!r}")
    return [int(v) for v in values]


def fractions(cell: str, field: str, number: int) -> list[Decimal]:
    values = parse_list(cell, field, number)
    require(all(type(v) in (int, float) and math.isfinite(v) and 0 <= v <= 1 for v in values),
            f"Data record {number}: invalid fractions in {field}: {cell!r}")
    require(math.isclose(math.fsum(values), 1.0, abs_tol=1e-6, rel_tol=0.0),
            f"Data record {number}: {field} sum differs from 1 beyond absolute tolerance 1e-6")
    return [Decimal(str(v)) for v in values]


def number(cell: str, field: str, row: int) -> Decimal:
    try:
        value = Decimal(cell.strip())
        require(value.is_finite() and math.isfinite(float(value)), "expected a finite numeric value")
        return value
    except (ValueError, InvalidOperation, OverflowError) as exc:
        raise ValueError(f"Data record {row}: invalid {field}: {cell!r}: {exc}") from exc


def observations(table: Table, binary: bool) -> list[Observation]:
    require(table.header == COLUMNS, f"{relative(table.path)}: original eight columns/order differ")
    result = []
    for index, row in enumerate(table.rows):
        record = dict(zip(table.header, row)); position = index + 1
        require(all(cell.strip() for cell in row), f"Data record {position}: missing/empty required cell")
        drug = ids(record["cmp_ids_drug"], "cmp_ids_drug", position)
        require(len(drug) == 1 and drug[0] in EXPECTED_COUNTS, f"Data record {position}: unexpected drug IDs {drug}")
        solvents = ids(record["cmp_ids_solvent"], "cmp_ids_solvent", position)
        composition = fractions(record["cmp_mole_fractions_solvent"], "cmp_mole_fractions_solvent", position)
        require(len(solvents) == len(composition) and len(solvents) == len(set(solvents)),
                f"Data record {position}: solvent IDs/fractions are not distinctly paired")
        active = sum(v > 0 for v in composition)
        require((len(solvents) == 2 and active == 2) if binary else active in (1, 2),
                f"Data record {position}: unexpected solvent composition, active={active}")
        drug_fractions = fractions(record["cmp_mole_fractions_drug"], "cmp_mole_fractions_drug", position)
        require(len(drug_fractions) == len(drug), f"Data record {position}: drug ID/fraction lengths differ")
        value, temperature = number(record["value"], "value", position), number(record["Temperature, K"], "Temperature, K", position)
        require(temperature > 0, f"Data record {position}: kelvin temperature is not greater than zero")
        result.append(Observation(index, row.copy(), drug[0], tuple(zip(solvents, composition)), value, temperature))
    return result


def constant(path: Path, name: str):
    tree = ast.parse(path.read_text(encoding="utf-8-sig"))
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == name for t in node.targets):
            return ast.literal_eval(node.value)
    raise ValueError(f"Constant {name} not found in {relative(path)}")


def find_unique(filename: str) -> Path:
    matches = [p for p in REPO_ROOT.rglob(filename) if p.is_file() and ".git" not in p.parts]
    require(len(matches) == 1, f"Expected one {filename}; found {[relative(p) for p in matches]}")
    return matches[0]


def citation(path: Path, text: str) -> str:
    for line, contents in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        if text in contents:
            return f"{relative(path)}:{line}"
    raise ValueError(f"Evidence text not found: {relative(path)}: {text}")


def ranges(positions: list[int]) -> str:
    groups = []
    for position in sorted(positions):
        if groups and groups[-1][1] + 1 == position:
            groups[-1][1] = position
        else:
            groups.append([position, position])
    return ", ".join(str(a) if a == b else f"{a}–{b}" for a, b in groups)


def duplicate_groups(models: list[Observation]) -> tuple[list[list[int]], list[list[int]], list[list[int]]]:
    exact, equivalent = defaultdict(list), defaultdict(list)
    for observation in models:
        exact[tuple(observation.row)].append(observation.index + 1)
        # Retain every other original cell literally, including value and drug lists.
        other_cells = tuple(cell for name, cell in zip(COLUMNS, observation.row)
                            if name not in {"cmp_ids_solvent", "cmp_mole_fractions_solvent"})
        equivalent[(other_cells, tuple(sorted(observation.pairs)))].append(observation.index + 1)
    exact_groups = [members for members in exact.values() if len(members) > 1]
    equivalent_groups = [members for members in equivalent.values() if len(members) > 1]
    order_only = [members for members in equivalent_groups
                  if len({tuple(cid for cid, _ in models[i - 1].pairs) for i in members}) > 1]
    return exact_groups, equivalent_groups, order_only


def validate_annotations(source: Table, models: list[Observation], annotated: Table) -> list[dict[str, str]]:
    require(annotated.header[:len(COLUMNS)] == COLUMNS, "Annotated original column order differs")
    required = {"source_row_index", "drug_id", "solvent_system_key", "reference_solvent_id",
                "reference_solvent_fraction", "other_solvent_id", "other_solvent_fraction", *FLAG_COLUMNS}
    require(required <= set(annotated.header), "Annotated traceability/grouping/flag columns are missing")
    require(len(annotated.rows) == len(source.rows), "Annotated rows added or missing")
    rows = annotated.dictionaries()
    require([row["source_row_index"] for row in rows] == list(map(str, range(len(source.rows)))),
            "source_row_index is not the unique zero-based input position in original order")
    for model, source_row, row in zip(models, source.rows, rows):
        require([row[name] for name in COLUMNS] == source_row, f"Annotated original cells changed at input record {model.index + 1}")
        require(row["drug_id"] == str(model.drug) and row["solvent_system_key"] == model.system,
                f"Annotated grouping differs at input record {model.index + 1}")
        pair = sorted(model.pairs)
        for prefix, (cid, fraction) in zip(("reference", "other"), pair):
            require(row[f"{prefix}_solvent_id"] == str(cid) and Decimal(row[f"{prefix}_solvent_fraction"]) == fraction,
                    f"Annotated ID/fraction alignment differs at input record {model.index + 1}")
    return rows


def reconcile_iqr(models: list[Observation], annotated: list[dict[str, str]], summary: Table, review: Table) -> dict:
    groups = defaultdict(list)
    for model in models:
        groups[("drug", str(model.drug), "")].append(model.index)
        groups[("drug_solvent_system", str(model.drug), model.system)].append(model.index)
    expected_keys = {(level, drug, system, variable) for level, drug, system in groups
                     for variable in ("value", "Temperature, K")}
    seen = set()
    with localcontext() as context:
        context.prec = 50
        for row in summary.dictionaries():
            key = (row["grouping_level"], row["drug_id"], row["solvent_system_key"], row["variable"])
            require(key in expected_keys and key not in seen, "IQR summary group/variable missing, unexpected or duplicated")
            seen.add(key)
            members = groups[key[:3]]
            require(int(row["observation_count"]) == len(members), f"IQR summary group count differs: {key}")
            variable = "value" if key[3] == "value" else "temperature"
            level = "drug" if key[0] == "drug" else "system"
            values = [models[i].value if variable == "value" else models[i].temperature for i in members]
            require(int(row["distinct_value_count"]) == len(set(values)) and Decimal(row["minimum"]) == min(values)
                    and Decimal(row["maximum"]) == max(values), f"IQR summary range/distinct count differs: {key}")
            eligible = level == "drug" or len(members) >= 10
            if not eligible:
                require(row["screening_status"] == "not screened—small group", "Small group incorrectly marked screened")
                require(all(annotated[i][f"{variable}_iqr_{d}_by_{level}"] == ""
                            for i in members for d in ("low", "high", "flag")), "Small-group flags are not missing")
                continue
            require(row["screening_status"] == "screened", "Eligible summary group is not screened")
            counts = {}
            for direction in ("low", "high", "flag"):
                column = f"{variable}_iqr_{direction}_by_{level}"
                require(all(annotated[i][column] in {"True", "False"} for i in members), "Invalid Boolean flag")
                counts[direction] = sum(annotated[i][column] == "True" for i in members)
            for i in members:
                low, high, flag = (annotated[i][f"{variable}_iqr_{d}_by_{level}"] == "True" for d in ("low", "high", "flag"))
                require(not (low and high) and flag == (low or high), "Directional flags disagree with combined flag")
            require((int(row["low_flag_count"]), int(row["high_flag_count"]), int(row["total_flagged_count"]))
                    == (counts["low"], counts["high"], counts["flag"]), f"IQR summary flag counts disagree: {key}")
            require(Decimal(row["flagged_percentage"]) == Decimal(counts["flag"]) * 100 / len(members), "IQR percentage denominator differs")
            # Check existing statistics algebra; do NOT calculate new quartiles/screens.
            q1, median, q3 = (Decimal(row[field]) for field in ("Q1", "median", "Q3"))
            iqr, lower, upper = (Decimal(row[field]) for field in ("IQR", "lower_fence", "upper_fence"))
            require(min(values) <= q1 <= median <= q3 <= max(values) and iqr == q3 - q1,
                    "Stored quartile/IQR algebra differs")
            require(lower == q1 - Decimal("1.5") * iqr and upper == q3 + Decimal("1.5") * iqr
                    and row["zero_iqr"] == str(iqr == 0), "Stored fence/zero-IQR algebra differs")
            for i, value in zip(members, values):
                if annotated[i][f"{variable}_iqr_low_by_{level}"] == "True":
                    require(value < lower, "Existing low flag does not cross its stored fence")
                if annotated[i][f"{variable}_iqr_high_by_{level}"] == "True":
                    require(value > upper, "Existing high flag does not cross its stored fence")
    require(seen == expected_keys, "IQR summary does not cover all comparison groups and variables")
    unions = {f"{variable}_by_{level}": {i for i, row in enumerate(annotated)
                                       if row[f"{variable}_iqr_flag_by_{level}"] == "True"}
              for level in ("drug", "system") for variable in ("value", "temperature")}
    union = set().union(*unions.values())
    review_rows = review.dictionaries()
    require(len(review_rows) == len(union) and {int(row["source_row_index"]) for row in review_rows} == union,
            "Review table is not the unique completed-screen union")
    for row in review_rows:
        index = int(row["source_row_index"])
        require(all(row[name] == annotated[index][name] for name in annotated[index]), "Review snapshot cells/flags differ")
        triggered = "; ".join(f"{variable}_by_{level}:{direction}"
                              for level in ("drug", "system") for variable in ("value", "temperature")
                              for direction in ("low", "high")
                              if row[f"{variable}_iqr_{direction}_by_{level}"] == "True")
        require(row["iqr_triggers"] == triggered, "Review trigger description differs")
    return {"counts": {key: len(indices) for key, indices in unions.items()},
            "flagged_indices": sorted(union), "review_count": len(review_rows), "summary_rows": len(seen)}


def verify_decision(input_path: Path, models: list[Observation], annotated: list[dict[str, str]],
                    decisions: Table, markdown: str, review: Table, source_dir: Path) -> dict[str, str]:
    rows = decisions.dictionaries()
    require(len(rows) == 1, "Expected exactly one authoritative decision")
    decision = rows[0]
    require(decision["source_row_index"] == "360" and len(models) > 360, "Decision identifies a different observation")
    model, flags = models[360], annotated[360]
    require(model.drug == 43 and dict(model.pairs) == {130: Decimal("0.945091"), 145: Decimal("0.054909")}
            and model.temperature == Decimal("273.15") and model.value == Decimal("-1.3493"), "Flagged observation identity differs")
    require(decision["input_dataset_path"] == relative(input_path) and decision["input_dataset_sha256"] == sha256(input_path),
            "Decision input path/hash does not identify the current binary input")
    require(decision["decision"] == "retain" and decision["review_status"] == "reviewed_retained", "Retention decision differs")
    require(decision["drug_id"] == "43" and decision["drug_name"] == "Posaconazole", "Decision drug identity differs")
    documented = {int(decision[f"solvent_id_{i}"]): Decimal(decision[f"solvent_fraction_{i}"]) for i in (1, 2)}
    require(documented == dict(model.pairs) and decision["solvent_system_key"] == model.system, "Decision solvent/fraction pairing differs")
    require(Decimal(decision["temperature_K"]) == model.temperature and Decimal(decision["stored_value"]) == model.value,
            "Decision measurement differs")
    require(decision["screening_level"] == "supplementary_by_system" and decision["variable"] == "value"
            and decision["flag_direction"] == "low", "Decision screening trigger differs")
    require(flags["value_iqr_low_by_system"] == flags["value_iqr_flag_by_system"] == "True"
            and flags["value_iqr_flag_by_drug"] == "False", "Retention has changed or cancelled statistical flags")
    for field in ("value_iqr_low_by_system", "value_iqr_flag_by_system", "value_iqr_flag_by_drug"):
        require(decision[field] == flags[field], "Decision recorded flags differ")
    require(decision["original_source_verification"] == "not independently verified during this review", "Source-verification limitation differs")
    for field in ("input_dataset_path", "input_dataset_sha256", "review_status", "decision", "review_method",
                  "decision_maker", "decision_reason", "original_source_verification", "supporting_figure_path",
                  "recorded_at", "recorded_timezone"):
        require(decision[field] in markdown, f"Decision Markdown disagrees with CSV: {field}")
    require("superseded by this separate decision record" in markdown, "Decision does not explain the screening snapshot")
    snapshot = [row for row in review.dictionaries() if row["source_row_index"] == "360"]
    require(len(snapshot) == 1 and snapshot[0]["review_status"] == "not reviewed", "Initial screening snapshot differs")
    require((REPO_ROOT / decision["supporting_figure_path"]).is_file(), "Supporting figure missing")
    names = {row["compound_id"]: row["name"].strip().casefold()
             for row in read_csv(source_dir / "processed_data/compounds.csv").dictionaries()}
    require(all(names.get(str(cid)) == name.casefold() for cid, name in {43: "Posaconazole", 130: "Ethanol", 145: "Ethyl acetate"}.items()),
            "Associated compound reference contradicts the reviewed identities")
    require(all(decision[f"solvent_name_{i}"].casefold() == names[decision[f"solvent_id_{i}"]]
                for i in (1, 2)), "Decision solvent display names disagree with explicit IDs")
    return decision


def historical_evidence(subset: Table, subset_models: list[Observation], binary_rows: list[list[str]], source_dir: Path) -> dict:
    source_name = constant(REPO_ROOT / "python scripts/05_subset_drugs_selection.py", "INPUT_NAME")
    full = read_csv(find_unique(source_name))
    require(full.header == COLUMNS, "Historical cumulative input schema differs")
    selected_indices = []
    for i, row in enumerate(full.rows):
        drug = ids(row[COLUMNS.index("cmp_ids_drug")], "cmp_ids_drug", i + 1)
        if len(drug) == 1 and drug[0] in EXPECTED_COUNTS:
            selected_indices.append(i)
    require([full.rows[i] for i in selected_indices] == subset.rows, "684-row subset does not equal the earlier complete ordered drug selection")
    mapping = read_csv(REPO_ROOT / "reports/data_cleaning/04_duplicate_removal/04_duplicate_row_mapping.csv").dictionaries()
    log = read_csv(REPO_ROOT / "reports/data_cleaning/03_repair/03_id_mismatch_change_log.csv").dictionaries()
    removal_log = read_csv(REPO_ROOT / "reports/data_cleaning/04_duplicate_removal/04_duplicate_removal_log.csv").dictionaries()
    repaired = read_csv(REPO_ROOT / "data/cleaned/drug-solubility/id_mismatch_repaired_DrugSolubilityData.csv")
    processed = read_csv(source_dir / "processed_data/processed_DrugSolubilityData.csv")
    raw = read_csv(source_dir / "raw_data/DrugSolubilityData.csv").dictionaries()
    compounds = {int(row["compound_id"]): row for row in read_csv(source_dir / "processed_data/compounds.csv").dictionaries()}
    smiles_ids = defaultdict(set)
    for cid, compound in compounds.items():
        smiles_ids[compound["smiles"]].add(cid)
    new_map = {int(row["new_index"]): row for row in mapping if row["status"] == "retained"}
    require(len(new_map) == len(full.rows) and len(repaired.rows) == len(processed.rows) == len(raw), "Historical row-map dimensions differ")
    logged = {int(row["data_row"]) - 1: row for row in log}
    changes, binary_changes = [], []
    for subset_index, full_index in enumerate(selected_indices):
        mapped = new_map[full_index]; original_index = int(mapped["original_index"])
        model = subset_models[subset_index]
        require(full.rows[full_index] == repaired.rows[original_index], "Historical retained row cells differ")
        require(raw[original_index]["Label"] == mapped["raw_label"], "Raw Label differs from the stable mapping")
        require(raw[original_index]["SMILES_0"] == compounds[model.drug]["smiles"], "Mapped raw drug identity differs")
        require(Decimal(raw[original_index]["LogS"]) == model.value
                and Decimal(raw[original_index]["Temperature, K"]) == model.temperature, "Mapped raw target/temperature differs")
        raw_pairs = []
        for slot in (1, 2):
            raw_fraction = Decimal(raw[original_index][f"comp_{slot}"]) / 100
            if raw_fraction:
                candidates = smiles_ids[raw[original_index][f"SMILES_{slot}"]]
                require(len(candidates) == 1, "Raw nonzero solvent identity is ambiguous")
                raw_pairs.append((next(iter(candidates)), raw_fraction))
        actual = sorted((cid, fraction) for cid, fraction in model.pairs if fraction > 0)
        expected = sorted(raw_pairs)
        require(len(actual) == len(expected) and all(cid == other and abs(frac - raw_frac) <= Decimal("1e-12")
                for (cid, frac), (other, raw_frac) in zip(actual, expected)), "Raw ID/fraction correspondence differs")
        differing = [name for j, name in enumerate(COLUMNS) if repaired.rows[original_index][j] != processed.rows[original_index][j]]
        if differing:
            require(differing == ["cmp_ids_solvent"] and original_index in logged, "Historical change is not a logged solvent-ID-only repair")
            item = logged[original_index]
            require(item["before"] == processed.rows[original_index][4] and item["after"] == repaired.rows[original_index][4]
                    and item["raw_label"] == mapped["raw_label"], "Historical repair log does not match the selected record")
            changes.append({"subset_data_row": subset_index + 1, "original_data_row": original_index + 1,
                            "raw_label": mapped["raw_label"], "drug_id": model.drug,
                            "active_solvents": model.active,
                            "positive_solvent_ids": [cid for cid, fraction in model.pairs if fraction > 0]})
            if model.row in binary_rows:
                binary_changes.append(subset_index + 1)
    removed = [row for row in mapping if row["status"] == "removed_duplicate"]
    require({row["original_data_row"] for row in removed} == {row["removed_original_data_row"] for row in removal_log},
            "Historical removal map and log disagree")
    target_removed = [row for row in removed if ids(repaired.rows[int(row["original_index"])][6], "cmp_ids_drug", int(row["original_data_row"]))
                      in ([5], [43], [86])]
    return {"selected_source": relative(full.path), "mapped_records": len(selected_indices),
            "raw_LogS_preserved": len(selected_indices), "subset_id_repairs": len(changes),
            "binary_id_repairs": len(binary_changes), "repair_details": changes,
            "three_drug_duplicate_removals": len(target_removed)}


def inspect_target_evidence(source_dir: Path) -> dict[str, str]:
    processing, metadata = source_dir / "processed_data/data_processing.py", source_dir / "processed_data/croissant.json"
    text = processing.read_text(encoding="utf-8-sig")
    tree = ast.parse(text)
    logarithms = [node for node in ast.walk(tree) if isinstance(node, ast.Call)
                  and isinstance(node.func, ast.Attribute) and node.func.attr in {"log", "log10", "log2", "exp"}]
    require(not logarithms, "Upstream processing contains a logarithmic/exponential call requiring new source review")
    description = json.loads(metadata.read_text(encoding="utf-8-sig"))["description"]
    evidence = {"rename": citation(processing, 'df = df.rename(columns={"LogS": "value"})'),
            "property": citation(processing, 'df["property"] = "Log solubility"'),
            "unit": citation(processing, 'df["unit"] = "g/L"'),
            "composition_percent": citation(processing, 'df[frac_col] = df[frac_col] / 100'),
            "raw_read": citation(processing, 'df = pd.read_csv(path)'),
            "metadata": citation(metadata, '"description": "Dataset of drug solubility'),
            "metadata_description": description,
            "raw_header": relative(source_dir / "raw_data/DrugSolubilityData.csv") + ":1"}
    for label, path, text in (("dataset_readme", source_dir.parent / "README.md", "**Drug Solubility:**"),
                              ("root_readme", source_dir.parent.parent / "README.md", "**Drug Solubility:**"),
                              ("unit_schema", source_dir.parent / "README.md", "`unit`:")):
        if path.is_file():
            evidence[label] = citation(path, text)
    return evidence


def validate(checks: Checks, input_path: Path, source_dir: Path) -> dict:
    context = {"input": relative(input_path)}
    filter_script, iqr_script = (REPO_ROOT / "python scripts" / name for name in ("06_remove_pure_solvents.py", "07_identify_iqr_outliers.py"))
    paths = {"binary": input_path, "annotated": input_path.with_name(input_path.stem + "_iqr_flagged.csv"),
             "summary": REPO_ROOT / "reports/iqr/iqr_summary.csv", "review": REPO_ROOT / "reports/iqr/iqr_review_rows.csv",
             "iqr_report": REPO_ROOT / "reports/iqr/iqr_report.md", "decisions": REPO_ROOT / "reports/iqr/iqr_review_decisions.csv",
             "decision_md": REPO_ROOT / "reports/iqr/iqr_review_decisions.md", "filter_script": filter_script, "iqr_script": iqr_script}
    missing = [relative(p) for p in paths.values() if not p.is_file()]
    checks.compare("required_files_present", [], missing)
    if missing:
        return context
    subset_name = constant(filter_script, "INPUT_NAME")
    require(constant(filter_script, "SUM_TOLERANCE") == 1e-6, "Filtering tolerance changed")
    subset_path = checks.attempt("original_subset_discovery", "one existing input named by the filtering script",
                                lambda: find_unique(subset_name), citation(filter_script, "INPUT_NAME ="))
    context["filter_evidence"] = [citation(filter_script, "INPUT_NAME ="), citation(filter_script, "def classify_record"),
                                  citation(filter_script, "active = sum(value > 0"), citation(filter_script, "def filter_records")]
    source = read_csv(input_path)
    checks.compare("required_columns_and_order", COLUMNS, source.header)
    if source.header != COLUMNS:
        return context
    missingness = {name: sum(not row[j].strip() for row in source.rows) for j, name in enumerate(COLUMNS)}
    context["missingness"] = missingness
    checks.compare("empty_required_cells", 0, sum(missingness.values()), show(missingness))
    models = checks.attempt("binary_structural_validity", "finite measurements; positive kelvin; one target drug; two positive paired solvents; valid drug fractions",
                            lambda: observations(source, True), "Safe literal parsing; absolute fraction-sum tolerance 1e-6, relative tolerance 0; no normalization.")
    checks.compare("analysis_record_count", 612, len(source.rows))
    if models is None:
        checks.add("dependent_reconciliations", "validated identities", "not evaluated", "SKIPPED", "Failed structural prerequisite")
        return context
    context["counts"] = dict(Counter(model.drug for model in models))
    checks.compare("analysis_drug_counts", EXPECTED_COUNTS, context["counts"])
    context["metadata"] = {}
    for field in ("property", "unit", "cmp_mole_fractions_drug"):
        distinct = sorted({row[COLUMNS.index(field)] for row in source.rows})
        context["metadata"][field] = distinct
        checks.compare("metadata_constant_" + field, 1, len(distinct), "Actual distinct stored cells: " + show(distinct) + "; consistency does not verify measurement definitions.")
    exact, equivalent, order_only = duplicate_groups(models)
    context["duplicates"] = {"exact_groups": exact, "equivalent_solvent_pair_groups": equivalent, "solvent_order_only_groups": order_only}
    checks.compare("exact_duplicate_groups_original_cells", 0, len(exact), "Keys include all eight original cells and value; audit columns excluded. Input data-record groups: " + show(exact))
    checks.add("equivalent_solvent_order_records", 0, len(order_only), "PASS" if not equivalent else "UNRESOLVED",
               "Only paired solvent IDs/fractions are canonicalized in memory, without tolerance; all other cells remain literal. "
               "Equivalent groups: " + show(equivalent) + "; ordering-only groups: " + show(order_only) + ". Candidates are not proof of copied experiments.")
    if subset_path is not None:
        subset = read_csv(subset_path)
        subset_models = checks.attempt("original_subset_structure", "valid one/two-active-solvent records under the existing filter rule",
                                       lambda: observations(subset, False))
        if subset_models is not None:
            excluded = [model for model in subset_models if model.active == 1]
            retained = [model for model in subset_models if model.active == 2]
            context["scope"] = {"input_path": relative(subset_path), "input_count": len(subset.rows),
                                "excluded_count": len(excluded), "retained_count": len(retained),
                                "original_counts": dict(Counter(m.drug for m in subset_models)),
                                "excluded_counts": dict(Counter(m.drug for m in excluded)),
                                "excluded_positions": {cid: ranges([m.index + 1 for m in excluded if m.drug == cid]) for cid in EXPECTED_COUNTS},
                                "excluded_solvents": dict(Counter(str(m.drug) + ":" + str(next(cid for cid, f in m.pairs if f > 0)) for m in excluded))}
            checks.compare("scope_original_record_count", 684, len(subset.rows))
            checks.compare("scope_original_drug_counts", EXPECTED_ORIGINAL, context["scope"]["original_counts"])
            checks.compare("pure_scope_exclusions", 72, len(excluded), "Every excluded record has exactly one positive solvent; all other records have exactly two.")
            checks.compare("pure_scope_exclusions_by_drug", EXPECTED_EXCLUDED, context["scope"]["excluded_counts"])
            checks.compare("scope_retained_records_equal_binary", True, [m.row for m in retained] == source.rows,
                           "Compare complete decoded fields in retained input order; no outlier criterion or subtraction-only inference.")
            history_paths = [source_dir / "raw_data/DrugSolubilityData.csv", source_dir / "processed_data/processed_DrugSolubilityData.csv",
                             REPO_ROOT / "data/cleaned/drug-solubility/id_mismatch_repaired_DrugSolubilityData.csv",
                             REPO_ROOT / "reports/data_cleaning/03_repair/03_id_mismatch_change_log.csv",
                             REPO_ROOT / "reports/data_cleaning/04_duplicate_removal/04_duplicate_row_mapping.csv",
                             REPO_ROOT / "reports/data_cleaning/04_duplicate_removal/04_duplicate_removal_log.csv"]
            if all(p.is_file() for p in history_paths):
                history = checks.attempt("subset_specific_history_and_raw_correspondence", "content-verified stable mapping and logs",
                                         lambda: historical_evidence(subset, subset_models, source.rows, source_dir),
                                         "Raw Labels are checked separately from positions. Raw fraction serialization comparison uses the existing 1e-12 tolerance; duplicate tests use no tolerance.")
                if history is not None:
                    context["history"] = history
                    checks.add("historical_solvent_id_repairs", "count supported by selected-record mappings/logs",
                               {"original_684_subset": history["subset_id_repairs"], "final_binary": history["binary_id_repairs"]}, "PASS",
                               "Only cmp_ids_solvent changed. Details: " + show(history["repair_details"]))
                    checks.add("historical_duplicate_removals_three_drugs", "count supported by removal mapping/log",
                               history["three_drug_duplicate_removals"], "PASS", "Full-dataset totals are not attributed to this subset.")
            else:
                checks.unresolved("subset_specific_history", "not established", "Missing historical evidence: " + show([relative(p) for p in history_paths if not p.is_file()]))
    annotated = checks.attempt("annotated_preservation_and_traceability", "all original cells/order; unique zero-based source index; aligned grouping fields",
                              lambda: validate_annotations(source, models, read_csv(paths["annotated"])))
    if annotated is not None:
        review = read_csv(paths["review"])
        iqr = checks.attempt("existing_iqr_summary_and_review_consistency", "group counts, flags, percentages and unique review union agree",
                             lambda: reconcile_iqr(models, annotated, read_csv(paths["summary"]), review),
                             "No new quartiles, flags or screening procedure executed; only existing results reconciled.")
        if iqr is not None:
            context["iqr"] = iqr
            for field, expected in {"value_by_drug": 0, "temperature_by_drug": 0, "value_by_system": 1, "temperature_by_system": 0}.items():
                checks.compare("existing_iqr_unique_flags_" + field, expected, iqr["counts"][field])
            checks.compare("existing_review_indices", [360], iqr["flagged_indices"])
        decision = checks.attempt("authoritative_retention_decision", "one matching reviewed_retained/retain decision with matching input hash and unchanged flag",
                                 lambda: verify_decision(input_path, models, annotated, read_csv(paths["decisions"]),
                                                         paths["decision_md"].read_text(encoding="utf-8"), review, source_dir))
        if decision is not None:
            context["decision"] = decision
            checks.compare("reviewed_retained_observations", 1, 1)
            checks.unresolved("retained_observation_original_accuracy", decision["original_source_verification"],
                              "Visual retention does not establish accuracy, a confirmed false positive or measurement error.")
        checks.compare("observations_removed_during_iqr", 0, len(source.rows) - len(annotated), "Full original-cell/order preservation is checked separately; all observations remain.")
        context["iqr_removed"] = len(source.rows) - len(annotated)
    report_text = paths["iqr_report"].read_text(encoding="utf-8")
    checks.compare("iqr_method_documentation", True, "exclude the overall middle observation" in report_text
                   and "n >= 10" in report_text and constant(iqr_script, "MIN_SYSTEM_N") == 10,
                   citation(iqr_script, "def quartiles") + "; " + citation(iqr_script, "lower, upper =") + "; reports/iqr/iqr_report.md: Method section. Median of halves and strict fences; per-drug primary and per-system supplementary; n>=10 is a project convention.")
    evidence = checks.attempt("local_target_processing_evidence", "inspect upstream processing and metadata",
                              lambda: inspect_target_evidence(source_dir))
    if evidence is not None:
        context["target_evidence"] = evidence
        checks.add("target_local_operation", "identify demonstrable local target operation", "LogS renamed to value; labels assigned; no logarithmic/exponential call in local upstream code", "PASS",
                   evidence["rename"] + "; " + evidence["property"] + "; " + evidence["unit"])
        for field, observed in (("target_logarithm_base", "not established"), ("target_concentration_basis", "not established"),
                                ("target_precompilation_conversion", "not established")):
            checks.unresolved(field, observed, evidence["raw_header"] + "; " + evidence["metadata"] + ". Metadata description: "
                              + evidence["metadata_description"] + ". Log solubility/g/L labels and an example g/100g solvent do not establish a per-observation formula or denominator.")
    return context


# Approved source-review reporting text, recorded 2026-10-14 (America/Toronto).
# These historical findings/versions are maintained here, not read from an output
# report or fetched from the web. Numerical findings still come from validate().
SOURCE_REVIEW_FINDINGS = """### Four documentation questions: current status after source review

| Question | Current status |
|---|---|
| 1. What logarithm base was used? | Established from the original authors' published code. |
| 2. What concentration basis was used? | Original target definition established; inconsistent inherited g/L metadata documented; metadata correction pending team discussion. |
| 3. What transformations occurred before value was created? | Main processing documented; earlier record-by-record literature extraction and preparation not independently audited. |
| 4. Why was the flagged observation retained? | Retention decision complete; independent verification against the original experimental publication not performed. |

### 1. Logarithm base

The original authors' external notebook, “Unit conversion,” `calculate_logs()`, explicitly applies `np.log10(df_updated['Solubility (g/100g)'])`. This establishes base 10 from the published code, rather than from the paper alone. [2]

### 2. Concentration basis and inherited unit discrepancy

Bao et al., **Results → Dataset overview** and **Methods → Data collection**, report LogS with S in g/100 g and describe standardizing literature measurements. [1] In `calculate_logs()`, drug mass is divided by the combined mass of drug and both solvents, then multiplied by 100. The denominator is total solution mass. [2]

The documented target is **LogS: the base-10 logarithm of the numerical solubility expressed in g/100 g of solution.** Thus `value` contains a logarithmic target, not an ordinary concentration measurement. [1, 2]

The provided Chemixhub code's `name_processing()` renames `LogS` to `value`, assigns `property = "Log solubility"`, and assigns `unit = "g/L"`. No corresponding numerical target conversion to g/L is implemented. Its division by 100 concerns composition percentages, not the target. [3]

**Metadata correction is pending discussion at the next group meeting.** This documentation update preserves all numerical values and inherited CSV labels, including `g/L`; it performs no conversion or metadata correction. The assigned label does not establish a g/L basis for the stored logarithmic target. [2, 3] The supplied `../chemixhub/datasets/drug-solubility/processed_data/croissant.json:51` also gives the generic example `g/100g solvent`; that wording differs from the total-solution denominator in the authors' code and is not used to override it. The dataset READMEs describe fields generically and do not resolve this discrepancy.

### 3. Documented transformation chain and provenance limit

The paper describes collection and standardization of literature solubility measurements. [1] The notebook's “Unit conversion” section implements solvent mole/mass-fraction interconversion in `calculate_fractions()`, followed by this target calculation in `calculate_logs()`: **standardized solubility input in mol/mol → component masses using composition and molecular weights → drug mass fraction → g/100 g of solution → base-10 logarithm, LogS**. [2] Chemixhub subsequently renames `LogS` to `value`. [3]

The raw header `../chemixhub/datasets/drug-solubility/raw_data/DrugSolubilityData.csv:1` names `LogS`. Existing content-verified project mappings establish numerical equality of raw `LogS` and stored `value` for all {mapped_records} selected records. This is internal preservation evidence, not independent verification of every source experiment or historical conversion.

The main calculation is documented; earlier record-by-record literature extraction and preparation of the standardized input have **not** been independently audited. Original experimental tables/SI and record-level preparation records would be needed for that audit. This is a provenance limitation, not evidence that the measurements are wrong."""

SOURCE_REVIEW_REFERENCES_AND_EVIDENCE = """## References and supporting sources

Sources were inspected on **2026-10-14 (America/Toronto)** during the recorded source review. Normal validation reruns do not access or reverify external sources. The published article was used; no local `1147_paper.pdf` was located. All four required sources were accessible. No source files were downloaded to the project.

**[1] Original research publication.** Bao, Z.; Tom, G.; Cheng, A.; Watchorn, J.; Aspuru-Guzik, A.; Allen, C. (2024). [Towards the prediction of drug solubility in binary solvent mixtures at various temperatures using machine learning](https://doi.org/10.1186/s13321-024-00911-3). *Journal of Cheminformatics*, 16, Article 117. [Published full text](https://link.springer.com/article/10.1186/s13321-024-00911-3). Relevant locations: **Results → Dataset overview**; **Methods → Data collection**; and the standalone **Data availability** heading, which links the authors' [Solubility_ML repository](https://github.com/Christine-Allen-Lab/Solubility_ML). Data availability is not a subsection of Supplementary Information. The explicit logarithm base is confirmed by code [2].

**[2] External supporting source code linked by the original publication.** Christine-Allen-Lab (n.d.). [Dataset_processing_overview.ipynb, Solubility_ML](https://github.com/Christine-Allen-Lab/Solubility_ML/blob/main/Codes_and_datasets/1_Codes_and_datasets_original/Exp_1_Raw_Dataset/Dataset_processing_overview.ipynb). Inspected [version e1b5c16712d14a398bf7229d37fb7775f31ea166](https://github.com/Christine-Allen-Lab/Solubility_ML/blob/e1b5c16712d14a398bf7229d37fb7775f31ea166/Codes_and_datasets/1_Codes_and_datasets_original/Exp_1_Raw_Dataset/Dataset_processing_overview.ipynb), read as notebook JSON in memory. Relevant location: **Unit conversion**, `calculate_fractions()` and `calculate_logs()` (notebook cells 38 and 39, zero-based; section heading in cell 37). This is externally accessed author code, not a notebook originally supplied in the Chemixhub dataset.

**[3] Provided dataset processing code / provenance evidence.** Chemixhub (n.d.). `data_processing.py`, supplied drug-solubility dataset. Dataset-repository-relative path: `datasets/drug-solubility/processed_data/data_processing.py`. Actual local location, relative to our project root: `../chemixhub/datasets/drug-solubility/processed_data/data_processing.py`, established from `08_final_validation.py`'s `DEFAULT_SOURCE` and recorded provenance. [Original repository URL](https://github.com/chemcognition-lab/chemixhub/blob/main/datasets/drug-solubility/processed_data/data_processing.py); [supplied checkout version 356c16b99e38b56a91b5219d33db410fc372d538](https://github.com/chemcognition-lab/chemixhub/blob/356c16b99e38b56a91b5219d33db410fc372d538/datasets/drug-solubility/processed_data/data_processing.py). The inspected local copy matches that checkout's tracked source after CRLF/LF normalization; its actual byte SHA-256 is `938f1a604477dcde28ddcc720bbe1b0cfdd03f619d28b182a30bf4d231aa7b17`. Relevant function: `name_processing()`, including composition-fraction handling, target renaming and label assignments. The local provided copy is the evidence for our processing history. The current online `main` version was not inspected or assumed identical. This is neither independent literature evidence nor code written by our project team.

**[4] External methodological reference for outlier treatment.** Frost, J. (n.d.). [Guidelines for Removing and Handling Outliers in Data](https://statisticsbyjim.com/basics/remove-outliers/). *Statistics By Jim*. Relevant location: **Guidelines for Dealing with Outliers**. Used only for the exclusion/retention principle; the project plot and saved review record supply the observation-specific assessment. No article publication date was established.

## Internal project evidence

These project-generated records support our counts, preservation checks and visual-review decision; external publications are not cited as evidence for those project results.

- Cleaning and selection: [03_repair_solvent_pairs.py](../../python%20scripts/03_repair_solvent_pairs.py), [04_remove_duplicates.py](../../python%20scripts/04_remove_duplicates.py), [05_subset_drugs_selection.py](../../python%20scripts/05_subset_drugs_selection.py), and [06_remove_pure_solvents.py](../../python%20scripts/06_remove_pure_solvents.py).
- Earlier preservation evidence: [ID-repair change log](../data_cleaning/03_repair/03_id_mismatch_change_log.csv), [ID-repair validation summary](../data_cleaning/03_repair/03_id_mismatch_repair_validation_summary.txt), [complete duplicate-removal row mapping](../data_cleaning/04_duplicate_removal/04_duplicate_row_mapping.csv), and [duplicate-removal log](../data_cleaning/04_duplicate_removal/04_duplicate_removal_log.csv).
- Screening implementation and historical results: [07_identify_iqr_outliers.py](../../python%20scripts/07_identify_iqr_outliers.py), [IQR summary](../iqr/iqr_summary.csv), [IQR report](../iqr/iqr_report.md), [unchanged screening review snapshot](../iqr/iqr_review_rows.csv), and [annotated dataset](../../data/processed/subset_drugs_binary_only_DrugSolubilityData_iqr_flagged.csv).
- Observation-specific review: [composition–value plot](../iqr/figures/drug_43_solvents_130_145_composition.png), [current decision CSV](../iqr/iqr_review_decisions.csv), and [current decision Markdown](../iqr/iqr_review_decisions.md).
- Validation: [08_final_validation.py](../../python%20scripts/08_final_validation.py) and [historical check results](final_validation_checks.csv)."""

OUTLIER_TREATMENT_NOTE = """Frost's “Guidelines for Dealing with Outliers” supports retaining observations unless exclusion has a specific justification; statistical extremeness alone is insufficient. [4] This is external methodological guidance, not evidence of this measurement's accuracy. The project-specific assessment comes from the user's saved visual review of the [composition–value plot](../iqr/figures/drug_43_solvents_130_145_composition.png) and the [separate retention decision](../iqr/iqr_review_decisions.md)."""


def render_source_review(context: dict) -> list[str]:
    """Combine the recorded source review with this run's verified context."""
    mapped = context.get("history", {}).get("raw_LogS_preserved", "not established")
    lines = SOURCE_REVIEW_FINDINGS.format(mapped_records=mapped).splitlines()
    lines += ["", "### 4. Retention of source_row_index 360", ""]
    if "decision" in context:
        decision = context["decision"]
        count = sum(context.get("counts", {}).values())
        lines += [f"The saved decision is **{decision['decision']}** for {decision['drug_name']} "
                  f"ID {decision['drug_id']} in ethanol–ethyl acetate, {decision['temperature_K']} K, "
                  f"stored value {decision['stored_value']}. `source_row_index = {decision['source_row_index']}` "
                  f"is the zero-based position in the **{count}-row binary-only input**, not the original "
                  f"full dataset. The supplementary per-system low-value flag remains **{decision['value_iqr_low_by_system']}**; "
                  f"the primary per-drug flag is **{decision['value_iqr_flag_by_drug']}**.", "",
                  "The user reported that the point appeared consistent with the surrounding composition trend at "
                  "the same temperature. The [decision CSV](../iqr/iqr_review_decisions.csv), "
                  "[decision explanation](../iqr/iqr_review_decisions.md), and "
                  "[supporting plot](../iqr/figures/drug_43_solvents_130_145_composition.png) are the project evidence. "
                  "Frost supports the treatment principle that extremeness alone is not sufficient justification "
                  "for removal. [4] Neither that article nor this visual decision independently verifies the original "
                  "measurement. Its accuracy is unconfirmed, and the flag has not been cancelled or established "
                  "as a false positive."]
    else:
        lines.append("This run did not establish the saved observation's correspondence and retention decision; "
                     "see the failed or skipped checks. The recorded source-review status does not override a failed check.")
    return lines


def render_summary(checks: Checks, context: dict, input_hash: str, protected_count: int) -> str:
    failures = [row for row in checks.rows if row["status"] == "FAIL"]
    outcome = "FINAL VALIDATION NOT PASSED" if failures else "STRUCTURAL / PRESERVATION CHECKS PASS; TARGET DEFINITION ESTABLISHED; INHERITED UNIT METADATA CORRECTION PENDING"
    lines = ["# Final validation and data-cleaning summary", "", f"**{outcome}**", "",
             "This is final checking and documentation of the existing strictly binary analysis population. "
             "No filtering, correction, outlier removal, screening rerun or new cleaned CSV was performed.", "",
             "Documentation amendment: sources reviewed on **2026-10-14 (America/Toronto)** established the target "
             "definition and updated the four documentation questions below. That documentation-only amendment "
             "did not rerun structural validation. The original **37 PASS, 0 FAIL, 4 UNRESOLVED** record remains "
             "unchanged in the validation-check CSV. Numerical values, inherited CSV labels, screening results "
             "and saved review decisions are unchanged.", "",
             "## Analysis population", "", f"- Analysis dataset: `{context['input']}`.",
             f"- Input SHA-256: `{input_hash}`.",
             "- Companion annotation: `data/processed/subset_drugs_binary_only_DrugSolubilityData_iqr_flagged.csv`.",
             "- Current review outcomes: `reports/iqr/iqr_review_decisions.csv` and `.md`.", "",
             "| Drug ID | Drug | Observed records |", "|---|---|---:|"]
    for cid, name in DRUG_NAMES.items():
        lines.append(f"| {cid} | {name} | {context.get('counts', {}).get(cid, 'not established')} |")
    lines += [f"| Total | | {sum(context['counts'].values()) if 'counts' in context else 'not established'} |", "",
              "## Integrity and duplicates", ""]
    if "missingness" in context:
        lines.append("Empty required cells by original column: " + json.dumps(context["missingness"], sort_keys=True) + ".")
    lines += ["", "Structural checks cover finite stored values and positive finite kelvin temperatures, exactly one "
              "target drug, exactly two distinct solvent IDs with aligned strictly positive finite fractions <=1, "
              "absolute fraction-sum tolerance 1e-6 (relative tolerance zero), and aligned valid drug fractions. "
              "No fraction normalization, imputation or unit change occurs. Each failure is recorded in the check CSV.", ""]
    if "duplicates" in context:
        result = context["duplicates"]
        lines += [f"Exact duplicate groups across the original eight decoded cells: **{len(result['exact_groups'])}**.",
                  f"Equivalent groups after in-memory paired solvent canonicalization: **{len(result['equivalent_solvent_pair_groups'])}**; "
                  f"groups differing in solvent order: **{len(result['solvent_order_only_groups'])}**.",
                  "Audit identifiers are excluded from all duplicate keys. All other cells, including stored value, "
                  "remain literal; ID/fraction numeric comparisons use exact equality, no tolerance or rounding. "
                  "Repeated conditions alone do not establish erroneous copies. No candidate is automatically removed."]
    lines += ["", "## Verified research-scope selection", ""]
    if "scope" in context:
        scope = context["scope"]
        lines += [f"Original subset: `{scope['input_path']}`, identified from INPUT_NAME in the filtering script.",
                  f"Actual reconciliation: **{scope['input_count']} original - {scope['excluded_count']} pure-solvent scope exclusions "
                  f"= {scope['retained_count']} binary observations**.", "",
                  "| Drug ID | Original | Pure excluded | Binary retained | Excluded original-subset data records |",
                  "|---|---:|---:|---:|---|"]
        for cid in EXPECTED_COUNTS:
            lines.append(f"| {cid} | {scope['original_counts'].get(cid, 0)} | {scope['excluded_counts'].get(cid, 0)} | "
                         f"{scope['original_counts'].get(cid, 0) - scope['excluded_counts'].get(cid, 0)} | {scope['excluded_positions'][cid]} |")
        lines += ["", "Those locations are **1-based data records in the original 684-row subset**, excluding its header; "
                  "they are calculated from actual compositions, not used as selection instructions. Every exclusion "
                  f"has exactly one positive solvent fraction, including {scope['excluded_solvents'].get('43:145', 0)} "
                  "pure ethyl-acetate Posaconazole records. All remaining records have two active solvents, and their complete "
                  "original cells match the binary input in order. These are scope exclusions, not outlier removals.",
                  "", "Excluded counts by drug ID:positive solvent ID: " + json.dumps(scope["excluded_solvents"], sort_keys=True) + "."]
        lines += ["", "Filtering evidence: " + "; ".join(f"`{ref}`" for ref in context["filter_evidence"]) + "."]
    else:
        lines.append("Scope reconciliation was not established; see failed/skipped checks.")
    lines += ["", "## Existing IQR results and completed review", "",
              "The existing procedure uses median-of-halves quartiles, excluding the overall middle observation "
              "for odd n, and strict 1.5 x IQR fences. Primary comparisons are per drug; supplementary "
              "comparisons are per drug and sorted solvent pair for n >= 10, a project convention, not a "
              "universal statistical minimum. Zero-IQR groups use the unchanged strict rule. See "
              "`reports/iqr/iqr_report.md` (Method section) and `python scripts/07_identify_iqr_outliers.py` "
              "(`quartiles` and `make_screen` sections). The screening script was not executed or changed.", ""]
    if "iqr" in context:
        iqr = context["iqr"]
        lines += ["Existing unique flag counts: " + json.dumps(iqr["counts"], sort_keys=True) + ".",
                  f"Review table unique observations: **{iqr['review_count']}**; source_row_index values: {iqr['flagged_indices']}.",
                  f"Summary group-variable rows reconciled: **{iqr['summary_rows']}**. Stored group counts, "
                  "percentages, flag columns, ranges and fence algebra agree. No new quartiles or flags were calculated."]
    if "decision" in context:
        decision = context["decision"]
        lines += ["", "The separate decision is authoritative: **reviewed_retained / retain**, made by the "
                  "project analyst (user) after visual inspection of the composition–value plot and temperature series.",
                  "", f"Observation: **source_row_index {decision['source_row_index']}** (zero-based position in the "
                  f"{sum(context.get('counts', {}).values())}-row binary input, not in the original full dataset), "
                  f"{decision['drug_name']} ID {decision['drug_id']}, ethanol ID {decision['solvent_id_1']} fraction "
                  f"{decision['solvent_fraction_1']}, ethyl acetate ID {decision['solvent_id_2']} fraction "
                  f"{decision['solvent_fraction_2']}, {decision['temperature_K']} K, stored value {decision['stored_value']}.",
                  "", "The recorded input SHA-256 matches the current binary file. The supplementary low-value "
                  "flag remains **True**, and the primary flag remains **False**. The unchanged automatically "
                  "generated review snapshot's initial 'not reviewed' entry is superseded by the separate "
                  "decision record. Retention does not cancel a statistical flag.", "",
                  "> " + decision["decision_reason"], "",
                  "Its accuracy was **not independently verified against the original experimental source**. "
                  "It is neither proven correct nor a confirmed error or false-positive flag.", "", OUTLIER_TREATMENT_NOTE]
    lines += ["", ("**Zero observations were removed during IQR screening/review**, confirmed by complete original-cell/order preservation."
                   if context.get("iqr_removed") == 0 else "The outlier-stage removal count was **not established** because preservation was not verified."), "",
              "## Subset-specific earlier history", ""]
    if "history" in context:
        history = context["history"]
        lines += [f"Content-verified mapping covers all **{history['mapped_records']}** selected observations, "
                  "using `reports/data_cleaning/04_duplicate_removal/04_duplicate_row_mapping.csv`, the cumulative "
                  "dataset selected by script 05, saved ID-repaired and upstream processed tables, raw Labels, "
                  "raw drug SMILES, temperatures and target values. Labels were read, not inferred from positions.",
                  "", f"- Logged solvent-ID-only corrections in the original 684-row subset: **{history['subset_id_repairs']}**.",
                  f"- Such corrections in the retained 612-row binary population: **{history['binary_id_repairs']}**.",
                  f"- Earlier exact-duplicate removals involving these three drug IDs: **{history['three_drug_duplicate_removals']}**.",
                  "- Correction membership and scope are established by the record details below; values, "
                  "temperatures and fractions were not repaired.",
                  "- These counts are derived from the selected records, repair change log and duplicate-removal "
                  "mapping/log. Full-dataset totals are not assigned to the three-drug subset. Evidence: "
                  "`reports/data_cleaning/03_repair/03_id_mismatch_change_log.csv` and "
                  "`reports/data_cleaning/04_duplicate_removal/04_duplicate_removal_log.csv`.",
                  "- Raw nonzero solvent ID/fraction correspondence uses the existing 1e-12 serialization "
                  "tolerance only for provenance matching, not duplicate detection or measurement similarity.", "",
                  "Historical correction details: " + json.dumps(history["repair_details"], sort_keys=True) + "."]
    else:
        lines.append("Subset-specific historical repair/removal counts: **not established** from available evidence; "
                     "no full-dataset totals have been substituted.")
    lines += ["", "## Target provenance versus measurement definitions", ""]
    lines += render_source_review(context)
    if "metadata" in context:
        lines += ["", "Actual distinct stored metadata: " + json.dumps(context["metadata"], sort_keys=True) + ". "
                  "Constant metadata and valid drug fractions establish consistency, not physical definitions."]
    lines += ["", "## Historical validation outcomes and safeguards", "",
              "Original saved result: **37 PASS, 0 FAIL, 4 UNRESOLVED**. The unchanged "
              "[validation-check CSV](final_validation_checks.csv) preserves the original check definitions, "
              "statuses and evidence; the later source review does not relabel them as a new all-PASS result.", "",
              "Protected inputs/companions at the original validation start: 61 files, including the validator itself. "
              "The original before/after hashes excluded only that validator's two output reports; protected contents "
              "were checked before and after publication.", "",
              f"Current validation status totals: {json.dumps(dict(Counter(row['status'] for row in checks.rows)), sort_keys=True)}. "
              "The same validation checks are rerun; their historical documentation questions are retained separately "
              "from the approved source-review conclusions. A rerun does not undo those later findings.", "",
              f"This run protects {protected_count} pre-existing inputs/companions, including the validator itself. "
              "Their hashes are checked before and after report publication; existing datasets, screening artifacts, "
              "decisions, scripts and figures are not modified by validation.", "",
              "Failed checks in this run:", ""]
    lines += [f"- **{row['check_name']}**: {row['observed_result']} — {row['notes']}" for row in failures] or ["- None."]
    lines += ["", "The four original documentation questions were `retained_observation_original_accuracy`, "
              "`target_logarithm_base`, `target_concentration_basis`, and `target_precompilation_conversion`. "
              "Their original `UNRESOLVED` rows remain intact. The current status table above records the later "
              "source-review findings: the base and original basis are established, the inherited unit label remains "
              "inconsistent pending team discussion, the main transformation is documented but earlier record-level "
              "preparation is unaudited, and retention is complete without original-source accuracy verification.", "",
              "Earlier IQR reports, plots and machine-readable validation records remain unchanged historical "
              "artifacts. Their original statements that definitions were unresolved are superseded **for this "
              "documentation review** by the sourced target definition here; their counts, flags and review history "
              "are preserved. Chemical-structure validation is outside this reporting task.", "",
              SOURCE_REVIEW_REFERENCES_AND_EVIDENCE, "", "## Reproduction and outputs", "",
              "From the repository root:", "", "```powershell", COMMAND, "```", "",
              "**Reporting synchronization (2026-10-15, America/Toronto):** `render_summary()` now independently "
              "generates the approved source findings, references, limitations and current review outcome. "
              "Validation was rerun for this synchronization task using the unchanged checks. This fixed date "
              "records the synchronization event, not the date of every later execution.", "",
              "Normal reruns generate reports from the validation inputs plus source-review text maintained in "
              "this script; the existing output Markdown is never used as a report template. No live web access "
              "is required. `publish_reports()` reuses byte-identical reports and refuses any unexpected existing "
              "report difference. The approved Markdown was synchronized once after snapshot verification and "
              "candidate-diff review; no force-overwrite option was added.", "",
              "The provided local Chemixhub checkout is required at `../chemixhub/datasets/drug-solubility` "
              "by default, together with its dataset/root READMEs and the existing project preservation evidence. "
              "An existing matching checkout can be selected with `--source-dir`; no source data is downloaded.", "",
              "- `python scripts/08_final_validation.py`.",
              "- `reports/validation/final_validation_checks.csv` — unchanged machine-readable statuses, "
              "expectations, observations and evidence.",
              "- `reports/validation/final_cleaning_summary.md` — this generated document.", "",
              "The existing binary-only CSV remains the analysis dataset; the annotated CSV and separate decision "
              "files are companion audit records. No newly renamed cleaned CSV was created. No files were staged, "
              "committed or pushed; the existing branch and unrelated work were preserved."]
    return "\n".join(lines) + "\n"


def publish_reports(payloads: dict[Path, bytes], protected_hashes: dict[Path, str]) -> None:
    for path, payload in payloads.items():
        require(path.resolve().is_relative_to(REPO_ROOT.resolve()) and not path.is_symlink(), "Unsafe report destination")
        require(path not in protected_hashes, "Report destination aliases a protected input")
        require(not path.exists() or (path.is_file() and path.read_bytes() == payload),
                f"Existing report differs; inspect and preserve it rather than overwrite: {relative(path)}")
    require(all(p.is_file() and sha256(p) == digest for p, digest in protected_hashes.items()), "Protected files changed before publication")
    temporary, created = [], []
    try:
        for path, payload in payloads.items():
            path.parent.mkdir(parents=True, exist_ok=True)
            if path.exists():
                continue
            stage = path.parent / f".final-validation-{uuid.uuid4().hex}.tmp"
            with stage.open("xb") as handle:
                temporary.append(stage)
                handle.write(payload)
            require(stage.read_bytes() == payload, "Temporary report read-back differs")
            os.link(stage, path)
            created.append(path)
        require(all(path.read_bytes() == payload for path, payload in payloads.items()), "Saved report read-back differs")
        require(all(p.is_file() and sha256(p) == digest for p, digest in protected_hashes.items()), "Protected files changed after publication")
    except BaseException:
        for path in reversed(created):
            path.unlink(missing_ok=True)
        raise
    finally:
        for stage in temporary:
            stage.unlink(missing_ok=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--report-dir", type=Path, default=REPO_ROOT / "reports/validation")
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE)
    args = parser.parse_args(argv)
    try:
        state = confirm_repository()
        input_path, report_dir, source_dir = map(resolve, (args.input, args.report_dir, args.source_dir))
        require(input_path.name == DEFAULT_INPUT.name and input_path.is_relative_to(REPO_ROOT.resolve()), "Wrong analysis input filename/location")
        output_csv, output_md = report_dir / "final_validation_checks.csv", report_dir / "final_cleaning_summary.md"
        protected = {}
        for base in (REPO_ROOT, source_dir):
            for path in base.rglob("*"):
                if path.is_file() and ".git" not in path.parts and path.resolve() not in {output_csv.resolve(), output_md.resolve()}:
                    protected[path] = sha256(path)
        for path in (source_dir.parent / "README.md", source_dir.parent.parent / "README.md"):
            if path.is_file():
                protected[path] = sha256(path)
        checks = Checks()
        checks.add("repository_and_branch", REPOSITORY + " / " + BRANCH, state["repository"] + " / " + state["branch"], "PASS")
        context = {"input": relative(input_path)}
        try:
            context = validate(checks, input_path, source_dir)
        except (OSError, ValueError, KeyError, TypeError, InvalidOperation, SyntaxError) as exc:
            checks.add("validation_completed", "all applicable checks evaluated", str(exc).replace(str(REPO_ROOT), "."), "FAIL", "Subsequent checks were not established; no data correction performed.")
        unchanged = all(path.is_file() and sha256(path) == digest for path, digest in protected.items())
        checks.compare("protected_preexisting_contents_unchanged", True, unchanged, f"{len(protected)} inputs/companions; no existing file edits.")
        current_state = confirm_repository()
        checks.compare("git_head_branch_and_staging_unchanged", state, current_state)
        buffer = io.StringIO(newline="")
        writer = csv.DictWriter(buffer, fieldnames=CHECK_COLUMNS)
        writer.writeheader(); writer.writerows(checks.rows)
        markdown = render_summary(checks, context, sha256(input_path) if input_path.is_file() else "not available", len(protected))
        publish_reports({output_csv: buffer.getvalue().encode("utf-8"), output_md: markdown.encode("utf-8")}, protected)
        saved = read_csv(output_csv)
        require(saved.header == CHECK_COLUMNS and saved.dictionaries() == checks.rows, "Saved machine-readable checks differ")
        require(output_md.read_text(encoding="utf-8") == markdown, "Saved Markdown differs")
        require(confirm_repository() == state, "Git state changed during report publication")
        failed = [row for row in checks.rows if row["status"] == "FAIL"]
        print("FINAL VALIDATION NOT PASSED" if failed else "STRUCTURAL / PRESERVATION VALIDATION PASS; SOURCE REVIEW DOCUMENTED; UNIT METADATA CORRECTION PENDING")
        print("Validation check statuses (historical documentation questions retained): "
              + show(dict(Counter(row["status"] for row in checks.rows))))
        print("Final drug counts: " + show(context.get("counts", "not established")))
        print("Scope: " + show(context.get("scope", "not established")))
        print("Existing IQR results: " + show(context.get("iqr", "not established")))
        print("Subset-specific history: " + show(context.get("history", "not established")))
        for row in checks.rows:
            if row["status"] in {"FAIL", "UNRESOLVED"}:
                prefix = "Historical documentation check UNRESOLVED" if row["status"] == "UNRESOLVED" else "FAIL"
                print(prefix + ": " + row["check_name"] + ": " + row["observed_result"])
        print("Current source-review conclusions and remaining limitations are reported separately in the Markdown; no live source verification was repeated.")
        print("Saved and read back: " + relative(output_csv))
        print("Saved and read back: " + relative(output_md))
        print("Protected files unchanged; no observations changed/removed; no files staged, committed or pushed.")
        return 2 if failed else 0
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"Final validation stopped: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
