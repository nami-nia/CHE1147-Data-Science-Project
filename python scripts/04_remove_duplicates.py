#!/usr/bin/env python3
"""Our project's keep-first policy for exact eight-cell CSV rows.

Python 3.10+; standard library only. Run from the group repository root:

    py -B -m unittest discover -s tests -p "test_remove_duplicates.py" -v
    py -B "python scripts/04_remove_duplicates.py" --input-csv "data/cleaned/drug-solubility/id_mismatch_repaired_DrugSolubilityData.csv" --dataset-dir "../chemixhub/datasets/drug-solubility" --dry-run

After reviewing the dry run, supply ALL FOUR new destinations:

    py -B "python scripts/04_remove_duplicates.py" --input-csv "data/cleaned/drug-solubility/id_mismatch_repaired_DrugSolubilityData.csv" --dataset-dir "../chemixhub/datasets/drug-solubility" --output-csv "data/cleaned/drug-solubility/id_mismatch_repaired_duplicate_removed_DrugSolubilityData.csv" --row-mapping "reports/data_cleaning/04_duplicate_removal/04_duplicate_row_mapping.csv" --removal-log "reports/data_cleaning/04_duplicate_removal/04_duplicate_removal_log.csv" --validation-summary "reports/data_cleaning/04_duplicate_removal/04_duplicate_removal_validation_summary.txt"

Default: dry run, no files or directories written. Explicit --dry-run cannot
accompany save arguments. Existing destinations and protected inputs are never
overwritten. A save reads all three CSV outputs back before writing the summary.
If saved-file validation fails, the exclusively created outputs are left for
inspection, their status is reported, and the command exits unsuccessfully.

The input is the cumulative Step 3 CSV, never the original processed CSV.
Its hash and the change log must match the Step 3 summary. Preserved input
order supplies a candidate raw correspondence; every candidate is checked
against raw targets, temperatures, and compound/fraction identities. Raw Labels
are copied, not generated from positions. After deletion the mapping selects
raw records for the unchanged Step 2 checks; the full raw table is never
compared positionally with the shortened output.

Duplicate keys are tuples of the eight parsed CSV strings, INCLUDING value.
No sorting, rounding, component normalization, or conditions-only deletion is
performed. Numeric tolerances inside existing structural checks are unchanged
and are never used for duplicate detection. CSV quoting may be reserialized;
every retained cell string and relative record order must remain identical.

Default expectations verify the approved 27,166 -> 27,157 dataset and nine
Label pairs. --expected-* arguments support verified synthetic fixtures or an
explicitly reviewed different input; they only check results, never select rows.
Chemical structures, unit definitions, experimental origin, and split design
remain separate review tasks. No existing split is applied or remapped here.
Exit codes: 0 = all checks passed; 2 = input, policy, validation, or output error.
"""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import os
import re
import runpy
import subprocess
import sys
from collections import Counter, defaultdict
from contextlib import ExitStack
from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path
from types import FunctionType


REPO_ROOT = Path(__file__).resolve().parent.parent
DATASET_NAME = "id_mismatch_repaired_duplicate_removed_DrugSolubilityData.csv"
INPUT_NAME = "id_mismatch_repaired_DrugSolubilityData.csv"
DEFAULT_SOURCE = REPO_ROOT.parent / "chemixhub/datasets/drug-solubility"
STEP3_FOLDER = REPO_ROOT / "reports/data_cleaning/03_repair"
COLUMNS = [
    "value", "Temperature, K", "property", "unit", "cmp_ids_solvent",
    "cmp_mole_fractions_solvent", "cmp_ids_drug", "cmp_mole_fractions_drug",
]
MAPPING_FIELDS = [
    "original_data_row", "new_data_row", "raw_label", "status",
    "original_index", "new_index",
]
LOG_FIELDS = [
    "duplicate_group", "removed_original_data_row", "removed_original_index",
    "removed_raw_label", "retained_original_data_row", "retained_original_index",
    "retained_raw_label", "retained_new_data_row", "retained_new_index", "reason",
]
REASON = "Exact duplicate across the eight dataset columns; first occurrence retained."
TEST_COMMAND = 'py -B -m unittest discover -s tests -p "test_remove_duplicates.py" -v'


@dataclass
class Inputs:
    dataset_dir: Path
    input_csv: Path
    step3_summary: Path
    step3_log: Path
    validator: Path
    repair_script: Path

    @property
    def raw(self) -> Path:
        return self.dataset_dir / "raw_data/DrugSolubilityData.csv"

    @property
    def compounds(self) -> Path:
        return self.dataset_dir / "processed_data/compounds.csv"

    @property
    def metadata(self) -> Path:
        return self.dataset_dir / "processed_data/croissant.json"


@dataclass
class Outputs:
    dataset: Path
    mapping: Path
    log: Path
    summary: Path

    def paths(self) -> list[Path]:
        return [self.dataset, self.mapping, self.log, self.summary]


@dataclass
class Result:
    rows: list[dict[str, str]]
    mapping: list[dict[str, str]]
    log: list[dict[str, str]]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def relative(path: Path) -> str:
    return Path(os.path.relpath(path.resolve(), REPO_ROOT)).as_posix()


def sha256(path: Path) -> str:
    with path.open("rb") as handle:
        digest = hashlib.sha256()
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path, *, fields: list[str] | None = None,
             required: set[str] | None = None) -> list[dict[str, str]]:
    """Reject malformed quoting, duplicate headers, widths, and schema changes."""
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, strict=True)
        header = next(reader, [])
        require(bool(header) and len(header) == len(set(header)),
                f"{path.name}: empty or duplicate header")
        require(fields is None or header == fields,
                f"{path.name}: expected columns in this order: {fields}")
        require(not required or required <= set(header),
                f"{path.name}: missing required columns")
        rows = []
        for record_n, cells in enumerate(reader, 1):
            require(len(cells) == len(header),
                    f"{path.name}: incorrect width at data row {record_n}")
            rows.append(dict(zip(header, cells)))
    return rows


def row_key(row: dict[str, str]) -> tuple[str, ...]:
    require(set(row) == set(COLUMNS) and all(isinstance(v, str) for v in row.values()),
            "Dataset records must contain exactly the eight string-valued columns")
    return tuple(row[column] for column in COLUMNS)


def deduplicate(rows: list[dict[str, str]], raw_labels: list[str]) -> Result:
    """Retain first occurrences in input order, without mutating any input."""
    require(len(rows) == len(raw_labels), "Unverified source correspondence: Label count")
    require(all(isinstance(label, str) and label.strip() for label in raw_labels)
            and len(set(raw_labels)) == len(raw_labels),
            "Raw Labels must be nonempty, unique source identifiers")
    output, mapping, log = [], [], []
    first, groups = {}, {}
    for old_index, (row, label) in enumerate(zip(rows, raw_labels)):
        key = row_key(row)
        kept = key not in first
        new_index = len(output) if kept else None
        mapping.append({
            "original_data_row": str(old_index + 1),
            "new_data_row": str(new_index + 1) if kept else "",
            "raw_label": label, "status": "retained" if kept else "removed_duplicate",
            "original_index": str(old_index), "new_index": str(new_index) if kept else "",
        })
        if kept:
            first[key] = (old_index, new_index)
            output.append(row.copy())
        else:
            if key not in groups:
                groups[key] = f"exact_duplicate_{len(groups) + 1:03d}"
            old_first, new_first = first[key]
            log.append({
                "duplicate_group": groups[key],
                "removed_original_data_row": str(old_index + 1),
                "removed_original_index": str(old_index), "removed_raw_label": label,
                "retained_original_data_row": str(old_first + 1),
                "retained_original_index": str(old_first),
                "retained_raw_label": raw_labels[old_first],
                "retained_new_data_row": str(new_first + 1),
                "retained_new_index": str(new_first), "reason": REASON,
            })
    return Result(output, mapping, log)


def summary_output_hashes(summary_path: Path) -> dict[Path, str]:
    text = summary_path.read_text(encoding="utf-8-sig")
    section = re.search(r"^OUTPUT SHA256\r?\n(.*?)(?:\r?\n\r?\n|\Z)", text, re.M | re.S)
    require(section is not None, "Step 3 summary has no OUTPUT SHA256 evidence")
    hashes = {}
    for line in section.group(1).splitlines():
        match = re.fullmatch(r"(.+): ([a-f0-9]{64})", line)
        require(match is not None, "Malformed Step 3 output hash evidence")
        path = (REPO_ROOT / match.group(1)).resolve()
        require(path not in hashes, "Ambiguous Step 3 output hash evidence")
        hashes[path] = match.group(2)
    return hashes


def verify_sources(paths: Inputs, rows: list[dict[str, str]], repair: dict) -> list[dict]:
    """Check Step 3 provenance and EVERY preserved raw/input correspondence."""
    original_processed = paths.dataset_dir / "processed_data/processed_DrugSolubilityData.csv"
    require(paths.input_csv.resolve() != original_processed.resolve(),
            "The original processed CSV cannot be the input for duplicate removal")
    evidence = summary_output_hashes(paths.step3_summary)
    for path in (paths.input_csv, paths.step3_log):
        require(path.resolve() in evidence and sha256(path) == evidence[path.resolve()],
                f"Step 3 output identity/hash mismatch: {relative(path)}")
    raw = read_csv(paths.raw, required=repair["RAW_REQUIRED"] | {"Label"})
    compounds = read_csv(paths.compounds, required={"compound_id", "smiles"})
    require(bool(rows), "Cleaning input has no records")
    labels = [r["Label"] for r in raw]
    require(all(label.strip() for label in labels) and len(set(labels)) == len(labels),
            "Raw Labels are missing or ambiguous")
    checked = repair["repair_rows"](raw, rows, compounds)
    require(not checked.changes and checked.rows == rows,
            "Input still requires solvent-ID repairs; duplicate removal will not repair it")
    entries = read_csv(paths.step3_log, fields=repair["LOG_FIELDS"])
    seen, fingerprint_cache = set(), {}
    for entry in entries:
        row_n = int(entry["data_row"])
        require(1 <= row_n <= len(rows) and row_n not in seen,
                "Invalid or repeated Step 3 change-log position")
        seen.add(row_n)
        source, after = raw[row_n - 1], rows[row_n - 1]
        require(entry["raw_label"] == source["Label"]
                and entry["after"] == after["cmp_ids_solvent"],
                "Step 3 change log disagrees with source/input identity")
        before = after.copy()
        before["cmp_ids_solvent"] = entry["before"]
        reconstructed = repair["repair_rows"]([source], [before], compounds)
        require(len(reconstructed.changes) == 1 and reconstructed.rows == [after],
                "Step 3 logged repair cannot be independently reproduced in memory")
        wanted = {**reconstructed.changes[0], "data_row": row_n}
        require(all(entry[k] == str(v) for k, v in wanted.items()),
                "Step 3 change-log fields disagree with the preserved record")
        fingerprints = json.loads(entry["input_sha256"])
        require(isinstance(fingerprints, dict) and bool(fingerprints),
                "Step 3 change log lacks input fingerprints")
        for key, digest in fingerprints.items():
            source_path = (REPO_ROOT / key).resolve()
            if source_path not in fingerprint_cache:
                fingerprint_cache[source_path] = sha256(source_path)
            require(fingerprint_cache[source_path] == digest,
                    f"Step 3 protected-source hash differs: {relative(source_path)}")
    return raw


def validate_mapping(input_rows: list[dict[str, str]], raw: list[dict],
                     result: Result) -> list[dict]:
    """Audit coverage, order, unchanged cells, and first-occurrence links."""
    require(len(input_rows) == len(raw) == len(result.mapping),
            "Mapping must account for every verified input record")
    by_removed = {entry["removed_original_index"]: entry for entry in result.log}
    require(len(by_removed) == len(result.log), "Removal log repeats an occurrence")
    selected_raw, first, seen_groups = [], {}, {}
    removed = 0
    for old_index, (original, source, mapping) in enumerate(zip(input_rows, raw, result.mapping)):
        require(set(mapping) == set(MAPPING_FIELDS), "Invalid mapping schema")
        require(mapping["original_data_row"] == str(old_index + 1)
                and mapping["original_index"] == str(old_index)
                and mapping["raw_label"] == source["Label"],
                "Mapping changed an original position or raw Label")
        key = row_key(original)
        if mapping["status"] == "retained":
            new_index = len(selected_raw)
            require(key not in first, "A later exact occurrence was incorrectly retained")
            require(mapping["new_data_row"] == str(new_index + 1)
                    and mapping["new_index"] == str(new_index),
                    "Retained mapping destinations have a gap, repetition, or reordering")
            require(new_index < len(result.rows) and result.rows[new_index] == original,
                    "Retained output cells differ from the ID-repaired input")
            first[key] = (old_index, new_index)
            selected_raw.append(source.copy())
        else:
            require(mapping["status"] == "removed_duplicate" and key in first,
                    "A unique or first-occurrence record was incorrectly removed")
            require(mapping["new_data_row"] == mapping["new_index"] == "",
                    "Removed records must have empty new positions")
            entry = by_removed.get(str(old_index))
            require(entry is not None and set(entry) == set(LOG_FIELDS),
                    "A removed record lacks its unique complete log entry")
            old_first, new_first = first[key]
            if key not in seen_groups:
                seen_groups[key] = f"exact_duplicate_{len(seen_groups) + 1:03d}"
            wanted = [seen_groups[key], str(old_index + 1), str(old_index), source["Label"],
                      str(old_first + 1), str(old_first), raw[old_first]["Label"],
                      str(new_first + 1), str(new_first), REASON]
            require([entry[k] for k in LOG_FIELDS] == wanted,
                    "Removal log does not identify the actual first retained counterpart")
            require(row_key(result.rows[new_first]) == key,
                    "Removed/retained records differ across the eight columns")
            removed += 1
    require(len(selected_raw) == len(result.rows) and removed == len(result.log),
            "Dataset, mapping, and removal-log counts disagree")
    return selected_raw


def structural_checks(paths: Inputs, validator: dict, selected_raw: list[dict],
                      dataset_path: Path, memory_rows: list[dict] | None = None) -> list[dict]:
    """Reuse every Step 2 check, with raw records selected by verified mapping.

    For saved validation, only the RAW loader is adapted. The unchanged
    processed loader reads the saved CSV directly. No ERROR is suppressed.
    The old validator's raw-derived candidate counts are not output statistics.
    """
    original = validator["validate"]
    environment = original.__globals__.copy()
    disk_loader = environment["load_csv"]

    def aligned_loader(path: Path, required: set[str]) -> list[dict]:
        if path.resolve() == paths.raw.resolve():
            require(all(required <= set(row) for row in selected_raw), "Invalid mapped raw schema")
            return [row.copy() for row in selected_raw]
        if memory_rows is not None and path.resolve() == dataset_path.resolve():
            require(all(required <= set(row) for row in memory_rows), "Invalid in-memory schema")
            return [row.copy() for row in memory_rows]
        return disk_loader(path, required)

    environment["load_csv"] = aligned_loader
    check = FunctionType(original.__code__, environment, original.__name__,
                         original.__defaults__, original.__closure__)
    _, issues, _ = check(paths.raw, dataset_path, paths.compounds, paths.metadata, 1e-6)
    errors = [issue for issue in issues if issue["severity"] == "ERROR"]
    require(not errors, "Structural/source validation failed: " + repr(errors[:3]))
    return issues


def decimal_list(text: str) -> list[Decimal]:
    # Preserve decimal digits while parsing existing numeric list cells. This
    # supports Python literal list syntax without rounding through float first.
    tree = ast.parse(text, mode="eval").body
    require(isinstance(tree, ast.List), "Expected a numeric list for condition counts")
    try:
        return [Decimal(ast.get_source_segment(text, node).replace("_", ""))
                for node in tree.elts]
    except (InvalidOperation, AttributeError) as exc:
        raise ValueError("Unsupported numeric list in condition count") from exc


def output_statistics(rows: list[dict]) -> dict:
    exact = Counter(row_key(row) for row in rows)
    conditions = defaultdict(list)
    for row in rows:
        pairs = []
        for kind in ("drug", "solvent"):
            ids = decimal_list(row[f"cmp_ids_{kind}"])
            fractions = decimal_list(row[f"cmp_mole_fractions_{kind}"])
            require(len(ids) == len(fractions), "ID/fraction length mismatch in output")
            pairs.append(tuple(sorted((cid, fraction) for cid, fraction in zip(ids, fractions)
                                      if fraction != 0)))
        key = (*pairs, Decimal(row["Temperature, K"]), row["property"], row["unit"])
        conditions[key].append(Decimal(row["value"]))
    repeated = [values for values in conditions.values() if len(values) > 1]
    return {
        "records": len(rows), "columns": len(COLUMNS),
        "exact_duplicate_groups": sum(n > 1 for n in exact.values()),
        "additional_exact_occurrences": sum(n - 1 for n in exact.values()),
        "repeated_condition_groups": len(repeated),
        "repeated_condition_records": sum(map(len, repeated)),
        "same_solubility_condition_groups": sum(len(set(v)) == 1 for v in repeated),
        "different_solubility_condition_groups": sum(len(set(v)) > 1 for v in repeated),
    }


def check_expectations(rows: list[dict], result: Result, args: argparse.Namespace) -> None:
    observed = (len(rows), len({e["duplicate_group"] for e in result.log}),
                len(result.log), len(result.rows))
    expected = (args.expected_input_records, args.expected_duplicate_groups,
                args.expected_removed_records, args.expected_output_records)
    require(observed == expected, f"Approved count expectations differ: expected {expected}; found {observed}")
    for field, wanted in [("retained_raw_label", args.expected_retained_labels),
                          ("removed_raw_label", args.expected_removed_labels)]:
        require([e[field] for e in result.log] == wanted,
                f"Approved {field} expectations differ; no outputs may be saved")


def protected_snapshot(paths: Inputs) -> dict[Path, str]:
    protected = set()
    for folder in (REPO_ROOT, paths.dataset_dir):
        for path in folder.rglob("*"):
            if path.is_file() and ".git" not in path.relative_to(folder).parts:
                protected.add(path.resolve())
    protected.update(p.resolve() for p in (paths.input_csv, paths.raw, paths.compounds,
                                          paths.metadata, paths.step3_log, paths.step3_summary,
                                          paths.validator, paths.repair_script, Path(__file__)))
    return {path: sha256(path) for path in protected}


def check_protected(before: dict[Path, str]) -> None:
    changed = [relative(path) for path, digest in before.items()
               if not path.is_file() or sha256(path) != digest]
    require(not changed, f"Protected pre-existing files changed: {changed}")


def guard_outputs(outputs: Outputs, protected: dict[Path, str]) -> None:
    resolved = [p.resolve() for p in outputs.paths()]
    require(len(set(resolved)) == 4, "Four output destinations must be distinct")
    for path, target in zip(outputs.paths(), resolved):
        require(target not in protected, f"Output would overwrite protected input: {relative(path)}")
        require(not path.exists() and not path.is_symlink(),
                f"Output already exists; refusing to overwrite: {relative(path)}")


def write_csv(handle, fields: list[str], rows: list[dict]) -> None:
    writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="raise")
    writer.writeheader()
    writer.writerows(rows)


def validate_result(paths: Inputs, rows: list[dict], raw: list[dict], result: Result,
                    validator: dict, dataset_path: Path, *, saved: bool) -> tuple[dict, list]:
    selected = validate_mapping(rows, raw, result)
    issues = structural_checks(paths, validator, selected, dataset_path,
                               None if saved else result.rows)
    stats = output_statistics(result.rows)
    require(stats["additional_exact_occurrences"] == 0, "Exact duplicates remain in output")
    second = deduplicate(result.rows, [r["Label"] for r in selected])
    require(not second.log and second.rows == result.rows, "Second pass removes additional records")
    return {**stats, "mapping_records": len(result.mapping), "removed_records": len(result.log),
            "retained_mappings": len(selected), "removed_mappings": len(result.log),
            "structural_errors": 0, "solvent_length_mismatches": 0,
            "second_pass_removals": len(second.log)}, issues


def portable_command(arguments: list[str]) -> str:
    path_options = {"--input-csv", "--dataset-dir", "--step3-summary", "--step3-log",
                    "--validator", "--repair-script", "--output-csv", "--row-mapping",
                    "--removal-log", "--validation-summary"}
    portable, previous = [], None
    for token in arguments:
        portable.append(relative(Path(token)) if previous in path_options else token)
        previous = token
    return subprocess.list2cmdline(["py", "-B", "python scripts/04_remove_duplicates.py", *portable])


def make_summary(paths: Inputs, outputs: Outputs, before: dict[Path, str],
                 stats: dict, issues: list, result: Result, command: str) -> str:
    reviews = [i for i in issues if i["severity"] == "REVIEW" and i["check"] == "UNIT_SOURCE_REVIEW"]
    lines = [
        "SAVED ID-REPAIRED AND EXACT-DUPLICATE-REMOVED DATASET VALIDATION",
        "Status: ALL REQUIRED SAVED-FILE CHECKS PASS; unresolved source review remains",
        f"Run UTC: {datetime.now(timezone.utc).isoformat(timespec='seconds')}",
        "", "PROJECT DECISION",
        "Our project's keep-first exact-row policy compares parsed CSV strings across all eight columns, including value.",
        REASON,
        "First means the earliest occurrence in the ID-repaired input order; no sorting or experimental ranking.",
        "This policy does not establish that the original researchers measured an experiment only once.",
        "No rounding, similarity tolerance, component normalization, averaging, unit conversion, or further ID repair was applied.",
        "CSV quoting is reserialized; retained cell strings and their relative order are unchanged.",
        "", "INPUT IDENTITY AND MAPPING",
        f"Cleaning input: {relative(paths.input_csv)}",
        "Input and Step 3 change-log hashes verified against Step 3 OUTPUT SHA256 evidence.",
        "Every input record checked against its order-preserved raw counterpart for target, temperature, drug, and solvent/fraction identities.",
        "Raw Labels are verified source fields, distinct from both old and new row positions.",
        "", "SAVED-FILE CHECKS (CSV FILES READ BACK FROM DISK)",
        f"Input records: {len(result.mapping):,}",
        f"Exact duplicate groups in input: {len({e['duplicate_group'] for e in result.log})}",
    ]
    lines.extend(f"{key}: {value:,}" for key, value in stats.items())
    lines += [
        "Original eight columns and column order: PASS",
        "Retained cells identical to corresponding ID-repaired input cells: PASS",
        "Retained relative order and first-occurrence policy: PASS",
        "Mapping covers every input row/index and every output row/index exactly once: PASS",
        "Removed mappings have empty new_data_row and new_index: PASS",
        "Each removed occurrence exactly matches its first retained counterpart in the saved output: PASS",
        "All prior solvent-ID repairs in retained records preserved: PASS",
        "Unchanged Step 2 checks reused with mapping-selected raw rows; no errors suppressed.",
        "The processed loader read the saved dataset directly; only raw selection was adapted in memory.",
        "Exact duplicates and repeated conditions recalculated from saved output cells, not the full raw table's old candidate report.",
        "Condition keys pair IDs with fractions, ignore only exact zero fractions, and include temperature/property/unit; value is evaluated afterward.",
        "No shortened-output positional comparison with the full raw table was performed.",
        "", "REMOVED -> RETAINED RAW LABELS",
    ]
    lines.extend(f"{e['removed_raw_label']} -> {e['retained_raw_label']} "
                 f"(old rows {e['removed_original_data_row']} -> {e['retained_original_data_row']}; "
                 f"retained new row {e['retained_new_data_row']})" for e in result.log)
    lines += ["", "INPUT / IMPLEMENTATION SHA256"]
    for path in (paths.input_csv, paths.raw, paths.compounds, paths.metadata,
                 paths.step3_summary, paths.step3_log, paths.validator, paths.repair_script,
                 Path(__file__), REPO_ROOT / "tests/test_remove_duplicates.py"):
        if path.is_file():
            lines.append(f"{relative(path)}: {sha256(path)}")
    lines += ["", "OUTPUT SHA256"]
    lines.extend(f"{relative(path)}: {sha256(path)}" for path in outputs.paths()[:3])
    lines += ["The summary has no self-referential hash; its final hash is printed to the console.",
              "", "PROTECTED PRE-EXISTING FILES: BEFORE AND AFTER"]
    lines.extend(f"{relative(path)}: before={digest}; after={sha256(path)}; unchanged=YES"
                 for path, digest in sorted(before.items(), key=lambda item: relative(item[0])))
    lines += ["", "UNRESOLVED LIMITATIONS"]
    lines.extend(f"{i['check']}: {i['detail']}" for i in reviews)
    lines += [
        "Chemical-structure validation is outside this task; no structures or compound references changed.",
        "Experimental measurement origin remains unresolved; exact table equality is not proof of experimental reuse.",
        "Existing train/validation/test split files were not applied, modified, regenerated, or remapped.",
        "Their positional indices cannot be reused unchanged after deletion and require a separate reviewed decision.",
        "Raw prediction columns and train/test assignments are unchanged; tracking is in separate mapping/log files.",
        "The dataset is not described as fully cleaned.",
        "", "REPRODUCTION (FROM THE GROUP REPOSITORY ROOT)", TEST_COMMAND,
        portable_command(["--input-csv", relative(paths.input_csv), "--dataset-dir", relative(paths.dataset_dir),
                          "--step3-summary", relative(paths.step3_summary), "--step3-log", relative(paths.step3_log),
                          "--validator", relative(paths.validator), "--repair-script", relative(paths.repair_script), "--dry-run"]),
        "Exact save command used:", command,
        "Save commands refuse existing outputs; use four new destinations to repeat a save.",
        "", "OUTPUT LOCATIONS",
    ]
    lines.extend(relative(path) for path in outputs.paths())
    return "\n".join(lines) + "\n"


def save_outputs(paths: Inputs, outputs: Outputs, before: dict[Path, str],
                 input_rows: list[dict], raw: list[dict], result: Result,
                 validator: dict, args: argparse.Namespace, command: str) -> dict:
    guard_outputs(outputs, before)
    created = []
    try:
        with ExitStack() as stack:
            handles = []
            for path in outputs.paths():
                path.parent.mkdir(parents=True, exist_ok=True)
                handle = stack.enter_context(path.open("x", encoding="utf-8", newline=""))
                created.append(path)
                handles.append(handle)
            for handle, fields, rows in zip(handles[:3], [COLUMNS, MAPPING_FIELDS, LOG_FIELDS],
                                             [result.rows, result.mapping, result.log]):
                write_csv(handle, fields, rows)
                handle.close()
            saved = Result(read_csv(outputs.dataset, fields=COLUMNS),
                           read_csv(outputs.mapping, fields=MAPPING_FIELDS),
                           read_csv(outputs.log, fields=LOG_FIELDS))
            check_expectations(input_rows, saved, args)
            stats, issues = validate_result(paths, input_rows, raw, saved, validator,
                                           outputs.dataset, saved=True)
            check_protected(before)
            summary = make_summary(paths, outputs, before, stats, issues, saved, command)
            handles[3].write(summary)
        require(outputs.summary.read_text(encoding="utf-8") == summary,
                "Saved summary readback differs from validated evidence")
        check_protected(before)
        return stats
    except (OSError, ValueError, csv.Error) as exc:
        status = ", ".join(f"{relative(p)} ({p.stat().st_size} bytes)" for p in created if p.exists())
        raise ValueError(f"Saved-output operation failed: {exc}. Created outputs left for review: {status or 'none'}") from exc


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-csv", type=Path, default=REPO_ROOT / "data/cleaned/drug-solubility" / INPUT_NAME)
    p.add_argument("--dataset-dir", type=Path, default=DEFAULT_SOURCE)
    p.add_argument("--step3-summary", type=Path, default=STEP3_FOLDER / "03_id_mismatch_repair_validation_summary.txt")
    p.add_argument("--step3-log", type=Path, default=STEP3_FOLDER / "03_id_mismatch_change_log.csv")
    p.add_argument("--validator", type=Path, default=Path(__file__).with_name("02_validate_values.py"))
    p.add_argument("--repair-script", type=Path, default=Path(__file__).with_name("03_repair_solvent_pairs.py"))
    p.add_argument("--dry-run", action="store_true", help="No output files or directories; also the default")
    p.add_argument("--output-csv", type=Path, help=f"New cumulative CSV, normally {DATASET_NAME}")
    p.add_argument("--row-mapping", type=Path)
    p.add_argument("--removal-log", type=Path)
    p.add_argument("--validation-summary", type=Path)
    for name, default in [("input-records", 27166), ("duplicate-groups", 9),
                          ("removed-records", 9), ("output-records", 27157)]:
        p.add_argument(f"--expected-{name}", type=int, default=default)
    p.add_argument("--expected-retained-labels", nargs="*", default=[str(i) for i in range(7905, 7914)])
    p.add_argument("--expected-removed-labels", nargs="*", default=[str(i) for i in range(7995, 8004)])
    return p


def main(argv: list[str] | None = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    p = parser()
    args = p.parse_args(arguments)
    destinations = [args.output_csv, args.row_mapping, args.removal_log, args.validation_summary]
    if any(path is not None for path in destinations) and not all(path is not None for path in destinations):
        p.error("Supply --output-csv, --row-mapping, --removal-log and --validation-summary together")
    if args.dry_run and args.output_csv is not None:
        p.error("--dry-run cannot be combined with save arguments")
    paths = Inputs(args.dataset_dir, args.input_csv, args.step3_summary, args.step3_log,
                   args.validator, args.repair_script)
    try:
        before = protected_snapshot(paths)
        outputs = Outputs(*destinations) if args.output_csv is not None else None
        if outputs is not None:
            guard_outputs(outputs, before)
        repair = runpy.run_path(str(paths.repair_script))
        validator = runpy.run_path(str(paths.validator))
        rows = read_csv(paths.input_csv, fields=COLUMNS)
        raw = verify_sources(paths, rows, repair)
        structural_checks(paths, validator, raw, paths.input_csv)
        result = deduplicate(rows, [r["Label"] for r in raw])
        check_expectations(rows, result, args)
        stats, _ = validate_result(paths, rows, raw, result, validator, paths.input_csv, saved=False)
        check_protected(before)
        print("KEEP-FIRST EXACT-ROW POLICY")
        print("Mode: DRY RUN; no files or directories written" if outputs is None else "Mode: explicit separate outputs")
        print(f"Verified ID-repaired input records: {len(rows):,}; source correspondence: PASS")
        print("ID-repaired disk-input structural errors: 0")
        print("PROPOSED OUTPUT VALIDATION (IN MEMORY; no saved output read yet)")
        print(json.dumps(stats, indent=2))
        for entry in result.log:
            print(f"Remove raw Label {entry['removed_raw_label']} (old row {entry['removed_original_data_row']}); "
                  f"retain {entry['retained_raw_label']} (old row {entry['retained_original_data_row']}, "
                  f"new row {entry['retained_new_data_row']})")
        print("Expected counts and retained/removed Labels: PASS")
        print("Protected pre-existing file hashes: unchanged")
        print("Unit definitions, structures, experimental origin and split design remain unresolved; no split files used.")
        if outputs is not None:
            stats = save_outputs(paths, outputs, before, rows, raw, result, validator, args,
                                 portable_command(arguments))
            print("SAVED OUTPUT VALIDATION (dataset, mapping and log read back from disk): PASS")
            print(json.dumps(stats, indent=2))
            for path in outputs.paths():
                print(f"{relative(path)}: SHA256 {sha256(path)}")
        return 0
    except (OSError, ValueError, csv.Error, KeyError, SyntaxError) as exc:
        print(f"Duplicate removal stopped: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
