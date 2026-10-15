#!/usr/bin/env python3
"""Flag potential IQR outliers without removing or modifying observations.

Python 3.10+; analysis uses the standard library, including Decimal arithmetic.
Figures require Matplotlib 3.10.9 and NumPy 2.2.6. From the repository root:

    uv run --no-project --with matplotlib==3.10.9 --with numpy==2.2.6 python -B "python scripts/07_identify_iqr_outliers.py"

With those plotting dependencies already installed, python -B also works.
Relative --input, --output, --report-dir and --compounds arguments start at this
script's repository, independently of the current working directory. Defaults:
data/processed/subset_drugs_binary_only_DrugSolubilityData.csv,
data/processed/subset_drugs_binary_only_DrugSolubilityData_iqr_flagged.csv,
reports/iqr/, and the sibling drug-solubility compounds.csv for display names.

source_row_index is created from the ZERO-BASED position in the 612-row input,
before grouping. It is not an index or Label in the original full dataset.
Stored fields and record order are preserved. Identity/fraction cells are safely
parsed; solvent pairs are sorted only in analysis fields, keeping ID/fraction
alignment. Only stored value and temperature receive IQR screens.

Quartiles use medians of halves, excluding the middle observation for odd n.
Fences are Q1 - 1.5*IQR and Q3 + 1.5*IQR; equality is not flagged. Decimal
precision is derived from the input, with no rounding before comparisons.
Zero IQR uses exactly the same rule. Primary groups are separate drugs;
supplementary drug/solvent groups require n >= 10, a PROJECT convention.
Smaller groups retain missing supplementary flags and report their ranges.

All outputs are staged temporarily, verified, and published without overwriting.
Identical existing artifacts are reused; conflicting artifacts cause failure.
The final CSVs are read back and checked again. Matplotlib caches and temporary
fixtures stay in the temporary staging directory and are cleaned up. No Git
operations, cleaning, transformations, imputation or treatment decisions occur.
Exit codes: 0 = complete verified screening; 2 = validation/output failure.
"""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import math
import os
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation, localcontext
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
INPUT_NAME = "subset_drugs_binary_only_DrugSolubilityData.csv"
OUTPUT_NAME = "subset_drugs_binary_only_DrugSolubilityData_iqr_flagged.csv"
DEFAULT_INPUT = REPO_ROOT / "data/processed" / INPUT_NAME
DEFAULT_OUTPUT = REPO_ROOT / "data/processed" / OUTPUT_NAME
DEFAULT_REPORT_DIR = REPO_ROOT / "reports/iqr"
DEFAULT_COMPOUNDS = REPO_ROOT.parent / "chemixhub/datasets/drug-solubility/processed_data/compounds.csv"
EXPECTED_COUNTS = {5: 360, 43: 162, 86: 90}
DRUG_NAMES = {5: "Itraconazole", 43: "Posaconazole", 86: "Ketoconazole"}
VARIABLES = ("value", "Temperature, K")
MIN_SYSTEM_N = 10
SMALL_STATUS = "not screened—small group"
ORIGINAL_COLUMNS = ["value", "Temperature, K", "property", "unit", "cmp_ids_solvent",
                    "cmp_mole_fractions_solvent", "cmp_ids_drug", "cmp_mole_fractions_drug"]
FLAG_COLUMNS = [f"{variable}_iqr_{direction}_by_{level}"
                for level in ("drug", "system") for variable in ("value", "temperature")
                for direction in ("low", "high", "flag")]
ADDED_COLUMNS = ["source_row_index", "drug_id", "solvent_system_key", "reference_solvent_id",
                 "reference_solvent_fraction", "other_solvent_id", "other_solvent_fraction",
                 "system_iqr_screening_status", *FLAG_COLUMNS]
REVIEW_COLUMNS = [*ORIGINAL_COLUMNS, *ADDED_COLUMNS, "iqr_triggers", "review_status"]
SUMMARY_COLUMNS = ["grouping_level", "drug_id", "solvent_system_key", "variable",
                   "observation_count", "distinct_value_count", "minimum", "Q1", "median", "Q3",
                   "maximum", "IQR", "lower_fence", "upper_fence", "low_flag_count",
                   "high_flag_count", "total_flagged_count", "flagged_percentage", "zero_iqr",
                   "screening_status"]
RUN_COMMAND = ('uv run --no-project --with matplotlib==3.10.9 --with numpy==2.2.6 '
               'python -B "python scripts/07_identify_iqr_outliers.py"')


@dataclass
class Record:
    source_row_index: int
    original: list[str]
    drug_id: int
    solvent_pair: tuple[int, int]
    fractions: dict[int, Decimal]
    value: Decimal
    temperature: Decimal

    @property
    def system_key(self) -> str:
        return "|".join(map(str, self.solvent_pair))


@dataclass
class Screen:
    summary: dict[str, str]
    indices: list[int]
    low: set[int]
    high: set[int]
    screened: bool


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def resolve_path(path: Path) -> Path:
    return (path if path.is_absolute() else REPO_ROOT / path).resolve()


def relative_path(path: Path) -> str:
    return Path(os.path.relpath(path, REPO_ROOT)).as_posix()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path) -> tuple[list[str], list[list[str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, strict=True)
        header = next(reader, [])
        require(bool(header) and all(name.strip() for name in header), f"{path}: empty/unnamed header")
        require(len(header) == len(set(header)), f"{path}: duplicate column names")
        rows = []
        try:
            for number, row in enumerate(reader, 1):
                require(len(row) == len(header), f"{path}: data record {number} has incorrect field count")
                rows.append(row)
        except csv.Error as exc:
            raise ValueError(f"{path}: invalid quoting at data record {len(rows) + 1}: {exc}") from exc
    return header, rows


def parse_list(cell: str, column: str, number: int) -> list:
    try:
        values = ast.literal_eval(cell)
        require(type(values) is list and bool(values), "expected a nonempty list")
        return values
    except (ValueError, SyntaxError, TypeError, OverflowError, RecursionError, MemoryError) as exc:
        raise ValueError(f"Data record {number}: invalid {column} cell {cell!r}: {exc}") from exc


def parse_ids(cell: str, column: str, number: int) -> list[int]:
    values = parse_list(cell, column, number)
    require(all(type(v) in (int, float) and v >= 0
                and (type(v) is int or (math.isfinite(v) and v.is_integer())) for v in values),
            f"Data record {number}: unsupported IDs in {column}: {cell!r}")
    return [int(v) for v in values]


def parse_fractions(cell: str, column: str, number: int, strictly_positive: bool) -> list[Decimal]:
    values = parse_list(cell, column, number)
    require(all(type(v) in (int, float) and math.isfinite(v)
                and (0 < v <= 1 if strictly_positive else 0 <= v <= 1) for v in values),
            f"Data record {number}: invalid fractions in {column}: {cell!r}")
    require(math.isclose(math.fsum(values), 1.0, rel_tol=0.0, abs_tol=1e-6),
            f"Data record {number}: {column} does not sum to 1 within absolute tolerance 1e-6: {cell!r}")
    return [Decimal(str(v)) for v in values]


def parse_number(cell: str, column: str, number: int) -> Decimal:
    try:
        value = Decimal(cell.strip())
        require(value.is_finite() and math.isfinite(float(value)), "expected a finite numeric value")
        return value
    except (InvalidOperation, ValueError, OverflowError) as exc:
        raise ValueError(f"Data record {number}: invalid {column} cell {cell!r}: {exc}") from exc


def parse_records(header: list[str], rows: list[list[str]]) -> list[Record]:
    require(header == ORIGINAL_COLUMNS, f"Unexpected input schema/order: {header}")
    columns = {name: header.index(name) for name in header}
    records = []
    for index, row in enumerate(rows):
        number = index + 1
        require(len(row) == len(header), f"Data record {number}: incorrect field count")
        drug = parse_ids(row[columns["cmp_ids_drug"]], "cmp_ids_drug", number)
        require(len(drug) == 1 and drug[0] in EXPECTED_COUNTS,
                f"Data record {number}: expected exactly one drug ID in {{5,43,86}}: {drug}")
        ids = parse_ids(row[columns["cmp_ids_solvent"]], "cmp_ids_solvent", number)
        fractions = parse_fractions(row[columns["cmp_mole_fractions_solvent"]],
                                    "cmp_mole_fractions_solvent", number, True)
        require(len(ids) == len(fractions) == 2 and len(set(ids)) == 2,
                f"Data record {number}: expected two distinct, correctly paired solvents")
        drug_fractions = parse_fractions(row[columns["cmp_mole_fractions_drug"]],
                                         "cmp_mole_fractions_drug", number, False)
        require(len(drug_fractions) == len(drug), f"Data record {number}: drug ID/fraction lengths differ")
        value = parse_number(row[columns["value"]], "value", number)
        temperature = parse_number(row[columns["Temperature, K"]], "Temperature, K", number)
        require(temperature > 0, f"Data record {number}: kelvin temperature must be greater than zero")
        # Create traceability before grouping. Sort IDs, never the fraction list.
        records.append(Record(index, row.copy(), drug[0], tuple(sorted(ids)),
                              dict(zip(ids, fractions)), value, temperature))
    return records


def median_sorted(values: list[Decimal]) -> Decimal:
    require(bool(values), "Cannot calculate a median from an empty group")
    middle = len(values) // 2
    return values[middle] if len(values) % 2 else (values[middle - 1] + values[middle]) / 2


def quartiles(values: list[Decimal]) -> tuple[Decimal, Decimal, Decimal]:
    require(len(values) >= 2, "Median-of-halves quartiles require at least two observations")
    ordered = sorted(values)  # Analysis-only copy; source records never move.
    middle = len(ordered) // 2
    lower = ordered[:middle]
    upper = ordered[middle + len(ordered) % 2:]
    return median_sorted(lower), median_sorted(ordered), median_sorted(upper)


def numeric_precision(records: list[Record]) -> int:
    values = [v for record in records for v in (record.value, record.temperature)]
    return max(50, max(v.adjusted() for v in values) - min(v.as_tuple().exponent for v in values) + 16)


def make_screen(members: list[Record], variable: str, level: str) -> Screen:
    values = [r.value if variable == "value" else r.temperature for r in members]
    ordered = sorted(values)
    eligible = level == "drug" or len(members) >= MIN_SYSTEM_N
    summary = {name: "" for name in SUMMARY_COLUMNS}
    summary.update(grouping_level=level, drug_id=str(members[0].drug_id),
                   solvent_system_key=members[0].system_key if level == "drug_solvent_system" else "",
                   variable=variable, observation_count=str(len(members)),
                   distinct_value_count=str(len(set(values))), minimum=str(ordered[0]),
                   median=str(median_sorted(ordered)), maximum=str(ordered[-1]),
                   screening_status="screened" if eligible else SMALL_STATUS)
    low, high = set(), set()
    if eligible:
        q1, q2, q3 = quartiles(values)
        iqr = q3 - q1
        lower, upper = q1 - Decimal("1.5") * iqr, q3 + Decimal("1.5") * iqr
        low = {r.source_row_index for r, v in zip(members, values) if v < lower}
        high = {r.source_row_index for r, v in zip(members, values) if v > upper}
        summary.update(Q1=str(q1), median=str(q2), Q3=str(q3), IQR=str(iqr),
                       lower_fence=str(lower), upper_fence=str(upper),
                       low_flag_count=str(len(low)), high_flag_count=str(len(high)),
                       total_flagged_count=str(len(low | high)),
                       flagged_percentage=str(Decimal(len(low | high)) * 100 / len(members)),
                       zero_iqr=str(iqr == 0))
    return Screen(summary, [r.source_row_index for r in members], low, high, eligible)


def triggers(row: dict[str, str]) -> str:
    return "; ".join(f"{variable}_by_{level}:{direction}"
                     for level in ("drug", "system") for variable in ("value", "temperature")
                     for direction in ("low", "high")
                     if row[f"{variable}_iqr_{direction}_by_{level}"] == "True")


def analyze(records: list[Record]) -> tuple[list[Screen], list[dict[str, str]], list[dict[str, str]]]:
    drugs, systems = defaultdict(list), defaultdict(list)
    annotated = []
    for record in records:
        drugs[record.drug_id].append(record)
        systems[(record.drug_id, record.solvent_pair)].append(record)
        reference, other = record.solvent_pair
        row = dict(zip(ORIGINAL_COLUMNS, record.original))
        row.update(source_row_index=str(record.source_row_index), drug_id=str(record.drug_id),
                   solvent_system_key=record.system_key, reference_solvent_id=str(reference),
                   reference_solvent_fraction=str(record.fractions[reference]), other_solvent_id=str(other),
                   other_solvent_fraction=str(record.fractions[other]), system_iqr_screening_status="")
        row.update({name: "" for name in FLAG_COLUMNS})
        annotated.append(row)
    screens = [make_screen(members, variable, level)
               for level, groups in (("drug", drugs), ("drug_solvent_system", systems))
               for _, members in sorted(groups.items()) for variable in VARIABLES]
    for screen in screens:
        variable = "value" if screen.summary["variable"] == "value" else "temperature"
        level = "drug" if screen.summary["grouping_level"] == "drug" else "system"
        for index in screen.indices:
            if level == "system":
                annotated[index]["system_iqr_screening_status"] = screen.summary["screening_status"]
            if screen.screened:
                for direction, flagged in (("low", index in screen.low), ("high", index in screen.high),
                                           ("flag", index in screen.low | screen.high)):
                    annotated[index][f"{variable}_iqr_{direction}_by_{level}"] = str(flagged)
    review = [dict(row, iqr_triggers=triggers(row), review_status="not reviewed")
              for row in annotated if triggers(row)]
    return screens, annotated, review


def self_checks() -> int:
    cases = [([1, 2, 3, 4, 5, 6, 7, 8], ("2.5", "4.5", "6.5"), 0, False),
             ([1, 2, 3, 4, 5, 6, 7, 8, 9], ("2.5", "5", "7.5"), 0, False),
             ([7] * 10, ("7", "7", "7"), 0, True),
             ([0] * 8 + [10], ("0", "0", "0"), 1, True),
             (["-1.5", 0, 0, 0, 1, 1, 1, "2.5"], ("0", "0.5", "1"), 0, False)]
    for values, expected, flagged, zero in cases:
        decimals = [Decimal(str(v)) for v in values]
        require(quartiles(decimals) == tuple(Decimal(v) for v in expected), "Quartile calculation check failed")
        members = [Record(i, [], 5, (127, 134), {}, value, Decimal("298.15"))
                   for i, value in enumerate(decimals)]
        screen = make_screen(members, "value", "drug")
        require(len(screen.low | screen.high) == flagged and screen.summary["zero_iqr"] == str(zero),
                "Fence/zero-IQR calculation check failed")
    row = ["1.00", "298.150", "Log solubility", "g/L", "[134,127.0]", "[0.2,0.8]", "[5.0]", "[1.0]"]
    record = parse_records(ORIGINAL_COLUMNS, [row])[0]
    require(record.source_row_index == 0 and record.solvent_pair == (127, 134)
            and record.fractions[127] == Decimal("0.8") and record.original == row,
            "Order-independent solvent pairing check failed")
    members = [Record(i, row.copy(), 5, (127, 134), record.fractions, Decimal(i), Decimal("298.15"))
               for i in range(9)]
    screens, annotated, review = analyze(members)
    require(all(s.summary["screening_status"] == SMALL_STATUS for s in screens if not s.screened)
            and all(r["value_iqr_flag_by_system"] == "" and r["temperature_iqr_flag_by_system"] == ""
                    for r in annotated), "Small-group flags must remain missing")
    return len(cases) + 2


def write_csv(path: Path, header: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def verify_saved(annotated_path: Path, summary_path: Path, review_path: Path,
                 records: list[Record], screens: list[Screen], expected_annotations: list[dict[str, str]]) -> None:
    header, cells = read_csv(annotated_path)
    require(header == ORIGINAL_COLUMNS + ADDED_COLUMNS, "Annotated schema/order changed")
    require(len(cells) == len(records), "Annotated record count changed")
    saved = [dict(zip(header, row)) for row in cells]
    require(saved == expected_annotations, "Annotated read-back differs from expected fields/flags/order")
    for record, row in zip(records, cells):
        require(row[:len(ORIGINAL_COLUMNS)] == record.original, "Original cell values changed")
        require(row[header.index("source_row_index")] == str(record.source_row_index), "Traceability/order changed")
    require(Counter(int(row["drug_id"]) for row in saved) == Counter(EXPECTED_COUNTS), "Annotated drug counts changed")
    summary_header, summary_cells = read_csv(summary_path)
    require(summary_header == SUMMARY_COLUMNS, "Summary schema changed")
    summaries = [dict(zip(summary_header, row)) for row in summary_cells]
    require(summaries == [screen.summary for screen in screens], "Summary read-back changed")
    for screen, row in zip(screens, summaries):
        variable = "value" if row["variable"] == "value" else "temperature"
        level = "drug" if row["grouping_level"] == "drug" else "system"
        require(int(row["observation_count"]) == len(screen.indices), "Summary group count mismatch")
        for direction, indices in (("low", screen.low), ("high", screen.high), ("flag", screen.low | screen.high)):
            column = f"{variable}_iqr_{direction}_by_{level}"
            if screen.screened:
                observed = {i for i in screen.indices if saved[i][column] == "True"}
                require(observed == indices, f"Summary/annotation {column} disagreement")
            else:
                require(all(saved[i][column] == "" for i in screen.indices), "Unscreened flags are not missing")
        if screen.screened:
            require(int(row["low_flag_count"]) == len(screen.low)
                    and int(row["high_flag_count"]) == len(screen.high)
                    and int(row["total_flagged_count"]) == len(screen.low | screen.high), "Summary flag totals mismatch")
            require(Decimal(row["flagged_percentage"]) == Decimal(len(screen.low | screen.high)) * 100 / len(screen.indices),
                    "Summary percentage uses an incorrect denominator")
    review_header, review_cells = read_csv(review_path)
    require(review_header == REVIEW_COLUMNS, "Review schema changed")
    expected_review = [dict(row, iqr_triggers=triggers(row), review_status="not reviewed")
                       for row in saved if triggers(row)]
    review_rows = [dict(zip(review_header, row)) for row in review_cells]
    require(review_rows == expected_review, "Review table is not exactly the unique union of completed flags")
    require(len({row["source_row_index"] for row in review_rows}) == len(review_rows), "Review rows duplicated")


def compound_names(path: Path) -> dict[int, str]:
    if not path.is_file():
        return {}
    header, rows = read_csv(path)
    require("compound_id" in header and "name" in header, "Compound display metadata lacks explicit ID/name fields")
    names = {}
    for number, row in enumerate(rows, 1):
        cid = parse_ids("[" + row[header.index("compound_id")] + "]", "compound_id", number)[0]
        name = row[header.index("name")].strip()
        if name:
            require(cid not in names or names[cid].casefold() == name.casefold(), f"Conflicting names for compound ID {cid}")
            names[cid] = name
    return names


def make_figures(directory: Path, records: list[Record], screens: list[Screen],
                 annotated: list[dict[str, str]], names: dict[int, str]) -> tuple[list[Path], dict[str, str]]:
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    import numpy
    from matplotlib.colors import Normalize

    directory.mkdir(parents=True)
    figures = []
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False}):
        for screen in screens:
            if screen.summary["grouping_level"] != "drug":
                continue
            cid, variable = int(screen.summary["drug_id"]), screen.summary["variable"]
            members = [records[i] for i in screen.indices]
            values = [r.value if variable == "value" else r.temperature for r in members]
            inliers = [v for r, v in zip(members, values) if r.source_row_index not in screen.low | screen.high]
            stats = {"label": DRUG_NAMES[cid], "q1": float(screen.summary["Q1"]),
                     "med": float(screen.summary["median"]), "q3": float(screen.summary["Q3"]),
                     "whislo": float(min(inliers)), "whishi": float(max(inliers)),
                     "fliers": [float(v) for r, v in zip(members, values) if r.source_row_index in screen.low | screen.high]}
            fig, ax = plt.subplots(figsize=(6.6, 5.2))
            ax.bxp([stats], showfliers=True, patch_artist=True,
                   boxprops={"facecolor": "#d1e5f0"}, medianprops={"color": "#16324f"},
                   flierprops={"marker": "x", "markeredgecolor": "#b2182b"})
            ax.axhline(float(screen.summary["lower_fence"]), color="#666666", linestyle="--", label="Lower IQR fence")
            ax.axhline(float(screen.summary["upper_fence"]), color="#999999", linestyle=":", label="Upper IQR fence")
            ax.set_ylabel("Stored value" if variable == "value" else "Temperature, K")
            ax.set_title(f"{DRUG_NAMES[cid]} (ID {cid}) — per-drug IQR\nMedian of halves; n = {len(members)}")
            ax.legend(loc="best", fontsize=9)
            ax.grid(axis="y", alpha=0.2)
            path = directory / f"drug_{cid}_{'value' if variable == 'value' else 'temperature'}_boxplot.png"
            fig.savefig(path, dpi=180, bbox_inches="tight", metadata={"Software": "CHE1147 IQR screening"})
            plt.close(fig)
            figures.append(path)

        systems = defaultdict(list)
        for record in records:
            systems[(record.drug_id, record.solvent_pair)].append(record)
        for (cid, pair), members in sorted(systems.items()):
            reference, other = pair
            temperatures = sorted({r.temperature for r in members})
            norm = Normalize(float(min(temperatures)), float(max(temperatures)))
            cmap = plt.get_cmap("viridis")
            fig, ax = plt.subplots(figsize=(8.2, 5.6))
            for temperature in temperatures:
                subset = sorted((r for r in members if r.temperature == temperature), key=lambda r: r.fractions[reference])
                ax.plot([float(r.fractions[reference]) for r in subset], [float(r.value) for r in subset],
                        color=cmap(norm(float(temperature))), alpha=0.5, linewidth=0.8)
            points = ax.scatter([float(r.fractions[reference]) for r in members], [float(r.value) for r in members],
                                c=[float(r.temperature) for r in members], cmap=cmap, norm=norm, s=30,
                                edgecolors="white", linewidths=0.3, label="Observations (color = temperature)")
            primary = [r for r in members if annotated[r.source_row_index]["value_iqr_flag_by_drug"] == "True"]
            supplementary = [r for r in members if annotated[r.source_row_index]["value_iqr_flag_by_system"] == "True"]
            if primary:
                ax.scatter([float(r.fractions[reference]) for r in primary], [float(r.value) for r in primary],
                           marker="x", s=110, c="#b2182b", linewidths=2, label="Primary per-drug value flag")
            if supplementary:
                ax.scatter([float(r.fractions[reference]) for r in supplementary], [float(r.value) for r in supplementary],
                           marker="D", s=110, facecolors="none", edgecolors="#e66101", linewidths=1.5,
                           label="Supplementary per-system value flag")
            bar = fig.colorbar(points, ax=ax, ticks=[float(t) for t in temperatures])
            bar.set_label("Temperature, K")
            ax.set_xlabel(f"Mole fraction of {names.get(reference, 'solvent')} (reference ID {reference})")
            ax.set_ylabel("Stored value")
            ax.set_xlim(0, 1)
            ax.set_title(f"{DRUG_NAMES[cid]} (ID {cid})\n{names.get(reference, str(reference))} (ID {reference}) + "
                         f"{names.get(other, str(other))} (ID {other})")
            ax.text(0.01, 0.02, f"Primary value flags: {len(primary)}; supplementary value flags: {len(supplementary)}",
                    transform=ax.transAxes, fontsize=8, bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"})
            ax.legend(loc="best", fontsize=8)
            ax.grid(alpha=0.15)
            path = directory / f"drug_{cid}_solvents_{reference}_{other}_composition.png"
            fig.savefig(path, dpi=180, bbox_inches="tight", metadata={"Software": "CHE1147 IQR screening"})
            plt.close(fig)
            figures.append(path)
    return figures, {"matplotlib": matplotlib.__version__, "numpy": numpy.__version__, "python": sys.version.split()[0]}


def display(value: str) -> str:
    if not value:
        return "—"
    number = Decimal(value)
    return format(number, ".6f").rstrip("0").rstrip(".") if number else "0"


def primary_overlap(annotated: list[dict[str, str]], cid: int | None = None) -> dict[str, int]:
    rows = [r for r in annotated if cid is None or int(r["drug_id"]) == cid]
    target = {r["source_row_index"] for r in rows if r["value_iqr_flag_by_drug"] == "True"}
    temperature = {r["source_row_index"] for r in rows if r["temperature_iqr_flag_by_drug"] == "True"}
    return {"target_only": len(target - temperature), "temperature_only": len(temperature - target),
            "both": len(target & temperature), "either": len(target | temperature)}


def write_report(path: Path, input_path: Path, annotated_path: Path, report_dir: Path, records: list[Record],
                 screens: list[Screen], annotated: list[dict[str, str]], review: list[dict[str, str]],
                 metadata: dict[str, list[str]], compounds: Path, versions: dict[str, str],
                 figure_names: list[str], input_hash: str, output_hashes: dict[str, str], checks: int) -> None:
    lines = ["# IQR screening of the strictly binary three-drug dataset", "",
             "All 612 observations remain in input order. Flags identify statistical review candidates; "
             "they are not confirmed errors. No source-based treatment decisions, removals, transformations "
             "or Git operations were made.", "", "## Input and verification", "",
             f"- Input: `{relative_path(input_path)}`; SHA-256 `{input_hash}`.",
             f"- Annotated output: `{relative_path(annotated_path)}`.",
             "- Validated counts: ID 5 = 360, ID 43 = 162, ID 86 = 90; total = 612.",
             "- Every record has one valid target drug, two distinct paired solvent IDs, and two finite "
             "strictly positive fractions no greater than 1. Solvent fractions sum to 1 with absolute "
             "tolerance 1e-6 and relative tolerance zero. Drug fractions are valid and aligned.",
             "- Stored value and kelvin temperature are finite numeric values; every temperature is positive.",
             "- Complete saved-file read-back checks passed: all original eight cells and column order, "
             "612-row order, source indices 0–611, drug counts, summary/flag agreement, and unique review union.",
             "- Input, original 684-record subset, other datasets and existing processing scripts were "
             "protected by before/after hashes and remained unchanged.",
             "- `source_row_index` is the original zero-based position in this 612-row INPUT ONLY. "
             "It is not a row position or raw Label in the full dataset.", "", "## Method", "",
             "For each group and variable, sort a temporary numerical copy. Q2 is the median of all values. "
             "For even n, split into equal halves. For odd n, exclude the overall middle observation. "
             "Q1 and Q3 are the medians of the lower and upper halves. IQR = Q3 − Q1; fences are "
             "Q1 − 1.5 × IQR and Q3 + 1.5 × IQR. Flag only strict inequalities; fence equality is retained "
             "without a flag. No interpolated percentile defaults are used.", "",
             "Decimal arithmetic uses precision derived from the input (at least 50 digits), with no "
             "rounding in flag comparisons. Tables in this report round only for display; the summary "
             "CSV retains the calculated quartiles and fences. Percentages use that particular group's n.", "",
             "If IQR is zero, the interval has zero width and the same strict rule applies, without "
             "epsilon. Constant groups have no flags. Nonconstant zero-IQR groups can flag departures "
             "from the common value and require caution.", "",
             "Primary comparisons combine all rows of each drug separately. Supplementary comparisons "
             "use drug plus an order-independent sorted solvent-ID pair. The n >= 10 requirement is "
             "a project screening convention, not a universal statistical minimum. Smaller groups "
             "are 'not screened—small group'; their supplementary flags are missing, not False. "
             "Their counts/ranges remain reported and their primary screens remain active.", "",
             "IDs are categorical. Composition is an experimental condition: no solvent-fraction IQR "
             "or automatic composition-based treatment was applied. The lower solvent ID is the "
             "reference for composition plots, with its own attached fraction.", "",
             f"Calculation checks: {checks} PASS (even, odd, constant, nonconstant zero-IQR, exact-fence, "
             "ID/fraction alignment and missing small-group flags).", "", "## Primary results", "",
             "| Drug ID | Variable | n | Distinct | Min | Q1 | Median | Q3 | Max | IQR | Lower fence | Upper fence | Low | High | Total | % | Zero IQR |",
             "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|"]
    fields = ["drug_id", "variable", "observation_count", "distinct_value_count", "minimum", "Q1", "median", "Q3",
              "maximum", "IQR", "lower_fence", "upper_fence", "low_flag_count", "high_flag_count", "total_flagged_count"]
    for screen in screens:
        if screen.summary["grouping_level"] == "drug":
            row = screen.summary
            lines.append("| " + " | ".join([row[f] if f in ("drug_id", "variable") else display(row[f]) for f in fields]
                                           + [f"{Decimal(row['flagged_percentage']):.2f}", row["zero_iqr"]]) + " |")
    lines += ["", "Temperature flags mean uncommon experimental temperatures within the comparison "
              "group; they are not automatically incorrect measurements or removal reasons.", "",
              "| Drug ID | Distinct temperatures, K |", "|---|---|"]
    for cid in EXPECTED_COUNTS:
        temperatures = sorted({r.temperature for r in records if r.drug_id == cid})
        lines.append(f"| {cid} | " + ", ".join(str(t) for t in temperatures) + " |")
    lines += ["", "Unique primary flagged rows (overlaps counted once):", "",
              "| Drug ID | Target only | Temperature only | Both | Either |", "|---|---:|---:|---:|---:|"]
    for cid in (*EXPECTED_COUNTS, None):
        overlap = primary_overlap(annotated, cid)
        lines.append(f"| {cid if cid is not None else 'All'} | " + " | ".join(str(v) for v in overlap.values()) + " |")
    lines += ["", "## Supplementary drug–solvent-system results", "",
              "Primary and supplementary flags are separate columns and distinct comparisons.", "",
              "| Drug ID | Solvent pair | Variable | n | Range | Q1 | Q3 | Fences | Low | High | Total | % | Status |",
              "|---|---|---|---:|---|---:|---:|---|---:|---:|---:|---:|---|"]
    for screen in screens:
        row = screen.summary
        if row["grouping_level"] == "drug_solvent_system":
            lines.append(f"| {row['drug_id']} | {row['solvent_system_key'].replace('|', '/')} | {row['variable']} | {row['observation_count']} | "
                         f"{display(row['minimum'])} to {display(row['maximum'])} | {display(row['Q1'])} | {display(row['Q3'])} | "
                         f"{display(row['lower_fence'])} to {display(row['upper_fence'])} | {row['low_flag_count']} | "
                         f"{row['high_flag_count']} | {row['total_flagged_count']} | {display(row['flagged_percentage'])} | {row['screening_status']} |")
    supplementary = {r["source_row_index"] for r in annotated
                     if r["value_iqr_flag_by_system"] == "True" or r["temperature_iqr_flag_by_system"] == "True"}
    lines += ["", f"Supplementary unique flagged rows: **{len(supplementary)}**. Unique union of ALL completed "
              f"screens in the review table: **{len(review)}**.", "",
              "| source_row_index | Input data record (1-based) | Drug ID | Solvent pair | Stored value | Temperature, K | Triggers | Review status |",
              "|---:|---:|---|---|---|---|---|---|"]
    for row in review:
        lines.append(f"| {row['source_row_index']} | {int(row['source_row_index']) + 1} | {row['drug_id']} | "
                     f"{row['solvent_system_key'].replace('|', '/')} | {row['value']} | {row['Temperature, K']} | "
                     f"{row['iqr_triggers']} | {row['review_status']} |")
    if not review:
        lines.append("| — | — | — | — | — | — | No completed screen flagged a row | — |")
    lines += ["", "## Metadata and unresolved measurement definitions", "",
              "| Field | Distinct stored cells | Constant? |", "|---|---|---|"]
    for field, values in metadata.items():
        lines.append(f"| {field} | `{json.dumps(values)}` | {len(values) == 1} |")
    lines += ["", "Drug-fraction cells were safely parsed and checked against the singleton drug lists; "
              "the confirmed constant field was not subjected to IQR. Neither metadata nor categorical IDs "
              "receive IQR flags.", "",
              "The local upstream `../chemixhub/datasets/drug-solubility/processed_data/data_processing.py` "
              "renames `LogS` to `value` and assigns `property = Log solubility` and `unit = g/L`. "
              "These operations and constant labels do not document the logarithm base or concentration "
              "denominator. Those definitions remain unresolved. Values were screened exactly as stored: "
              "no extra logarithm or back-transformation. Figures say **Stored value**.", "",
              f"Compound display-name source: `{relative_path(compounds)}` "
              f"({'available; explicit compound IDs' if compounds.is_file() else 'unavailable; IDs shown instead'}). "
              "Names are plot labels, not chemical-structure validation.", "",
              "Flags depend on comparison scope and can reflect experimental composition/temperature "
              "patterns. Inspect the composition plots and source evidence before any later treatment. "
              "No observation has been removed or declared erroneous; review_status remains not reviewed.", "",
              "## Figures and reproduction", "",
              "Boxplots use precomputed Q1, median and Q3 passed to Matplotlib bxp, not percentile defaults. "
              "Whiskers end at the most extreme observed values inside the exact fences; fence lines are "
              "drawn separately. Composition plots identify the lower-ID reference solvent, color all "
              "observations by temperature, mark primary value flags with red crosses and supplementary "
              "value flags with orange diamonds. Lines connect observed points at each temperature for "
              "context; they are not fitted models or treatment rules.", "",
              f"Dependencies used: Python {versions['python']}, Matplotlib {versions['matplotlib']}, NumPy {versions['numpy']}.",
              "Numerical analysis, validation and CSV handling otherwise use the standard library.", "",
              "Run from the repository root (defaults resolve relative to the script's repository):", "",
              "```powershell", RUN_COMMAND, "```", "",
              "Identical rerun artifacts are reused unchanged; any conflicting existing artifact is "
              "refused before publishing. Temporary staging and plotting caches are removed.", "",
              "Generated outputs:", "",
              f"- `{relative_path(annotated_path)}`",
              f"- `{relative_path(report_dir / 'iqr_summary.csv')}`",
              f"- `{relative_path(report_dir / 'iqr_review_rows.csv')}`",
              f"- `{relative_path(report_dir / 'iqr_report.md')}`"]
    lines += [f"- `{relative_path(report_dir / 'figures' / name)}`" for name in figure_names]
    lines += ["", "Saved CSV SHA-256 values:", ""]
    lines += [f"- `{name}`: `{digest}`" for name, digest in output_hashes.items()]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def guard_destinations(paths: list[Path], protected: list[Path]) -> None:
    require(len({p.resolve() for p in paths}) == len(paths), "Output destinations alias each other")
    for path in paths:
        require(not path.is_symlink() and path.resolve().is_relative_to(REPO_ROOT.resolve()),
                f"Output must be a regular path inside the repository: {path}")
        require(not path.exists() or path.is_file(), f"Output is not a regular file: {path}")
        for source in protected:
            require(path.resolve() != source.resolve(), f"Output aliases protected source: {path}")
            if path.exists() and source.exists():
                require(not os.path.samefile(path, source), f"Output aliases protected source: {path}")


def check_protected(hashes: dict[Path, str]) -> None:
    changed = [relative_path(p) for p, digest in hashes.items() if not p.is_file() or sha256(p) != digest]
    require(not changed, f"Protected source files changed: {changed}")


def publish_artifacts(artifacts: dict[Path, Path], protected: list[Path], hashes: dict[Path, str]) -> list[Path]:
    guard_destinations(list(artifacts), protected)
    # Check every collision before creating any final artifact.
    for destination, staged in artifacts.items():
        require(not destination.exists() or destination.read_bytes() == staged.read_bytes(),
                f"Existing output conflicts; refusing overwrite: {destination}")
    check_protected(hashes)
    created = []
    try:
        for destination, staged in artifacts.items():
            destination.parent.mkdir(parents=True, exist_ok=True)
            if destination.exists():
                require(destination.read_bytes() == staged.read_bytes(), f"Output changed before publication: {destination}")
            else:
                os.link(staged, destination)
                created.append(destination)
        check_protected(hashes)
        return created
    except BaseException:
        for path in reversed(created):
            path.unlink(missing_ok=True)
        raise


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--compounds", type=Path, default=DEFAULT_COMPOUNDS, help="Explicit compound-ID/name display metadata")
    args = parser.parse_args(argv)
    created = []
    try:
        input_path, output, report_dir, compounds = map(resolve_path, (args.input, args.output, args.report_dir, args.compounds))
        require(input_path.is_file() and input_path.name == INPUT_NAME, f"Required binary input missing/wrong filename: {input_path}")
        require(output.name == OUTPUT_NAME, f"Annotated output filename must be {OUTPUT_NAME}")
        require(input_path.is_relative_to(REPO_ROOT.resolve()), "Input must be in this repository")
        require(report_dir.is_relative_to(REPO_ROOT.resolve()), "Reports must stay inside the repository")
        protected = [input_path, Path(__file__).resolve(),
                     REPO_ROOT / "data/processed/subset_drugs_DrugSolubilityData.csv",
                     *sorted((REPO_ROOT / "data/cleaned").rglob("*.csv")),
                     *sorted((REPO_ROOT / "python scripts").glob("*.py"))]
        if compounds.is_file():
            protected.append(compounds)
        hashes = {p: sha256(p) for p in set(protected) if p.is_file()}
        header, rows = read_csv(input_path)
        records = parse_records(header, rows)
        counts = Counter(r.drug_id for r in records)
        print(f"Input: {relative_path(input_path)}; records={len(records)}; per-drug counts={dict(counts)}")
        require(counts == Counter(EXPECTED_COUNTS) and len(records) == 612,
                f"Input count discrepancy: actual={dict(counts)} total={len(records)}; expected={EXPECTED_COUNTS} total=612")
        metadata = {name: sorted({row[header.index(name)] for row in rows})
                    for name in ("property", "unit", "cmp_mole_fractions_drug")}
        names = compound_names(compounds)
        with localcontext() as context:
            context.prec = numeric_precision(records)
            checks = self_checks()
            print(f"Quartile, fence, pairing and small-group calculation checks: {checks} PASS")
            screens, annotated, review = analyze(records)
            with tempfile.TemporaryDirectory(prefix=".iqr-build-", dir=REPO_ROOT) as staging_name:
                staging = Path(staging_name)
                staged_annotated, staged_summary, staged_review = staging / OUTPUT_NAME, staging / "iqr_summary.csv", staging / "iqr_review_rows.csv"
                write_csv(staged_annotated, ORIGINAL_COLUMNS + ADDED_COLUMNS, annotated)
                write_csv(staged_summary, SUMMARY_COLUMNS, [screen.summary for screen in screens])
                write_csv(staged_review, REVIEW_COLUMNS, review)
                previous_config = os.environ.get("MPLCONFIGDIR")
                os.environ["MPLCONFIGDIR"] = str(staging / "matplotlib-cache")
                try:
                    figures, versions = make_figures(staging / "figures", records, screens, annotated, names)
                finally:
                    if previous_config is None:
                        os.environ.pop("MPLCONFIGDIR", None)
                    else:
                        os.environ["MPLCONFIGDIR"] = previous_config
                verify_saved(staged_annotated, staged_summary, staged_review, records, screens, annotated)
                check_protected(hashes)
                staged_report = staging / "iqr_report.md"
                output_hashes = {relative_path(output): sha256(staged_annotated),
                                 relative_path(report_dir / "iqr_summary.csv"): sha256(staged_summary),
                                 relative_path(report_dir / "iqr_review_rows.csv"): sha256(staged_review)}
                write_report(staged_report, input_path, output, report_dir, records, screens, annotated, review,
                             metadata, compounds, versions, [p.name for p in figures], hashes[input_path], output_hashes, checks)
                artifacts = {output: staged_annotated, report_dir / "iqr_summary.csv": staged_summary,
                             report_dir / "iqr_review_rows.csv": staged_review, report_dir / "iqr_report.md": staged_report}
                artifacts.update({report_dir / "figures" / p.name: p for p in figures})
                created = publish_artifacts(artifacts, protected, hashes)
                verify_saved(output, report_dir / "iqr_summary.csv", report_dir / "iqr_review_rows.csv", records, screens, annotated)
                require(all(destination.read_bytes() == staged.read_bytes() for destination, staged in artifacts.items()),
                        "Final artifact differs from its validated staged content")
                check_protected(hashes)
            print("COMPLETE SAVED-OUTPUT VERIFICATION: PASS")
            print("612 observations; original fields, columns and order preserved; source_row_index=0..611: PASS")
            print("Summary counts/percentages agree with flags; review is the unique completed-screen union: PASS")
            print("Input, original subset and protected datasets/scripts unchanged: PASS")
            for screen in screens:
                row = screen.summary
                if row["grouping_level"] == "drug" or row["total_flagged_count"] not in ("", "0"):
                    print(json.dumps(row))
            print("Unique primary overlaps: " + json.dumps(primary_overlap(annotated)))
            print(f"Unique review rows: {len(review)}; summary rows: {len(screens)}; figures: {len(figures)}")
            for path in artifacts:
                print(("Created: " if path in created else "Reused unchanged: ") + relative_path(path))
            print("No observations removed; no Git operations performed.")
            return 0
    except (OSError, ValueError, csv.Error, UnicodeError, InvalidOperation, ImportError) as exc:
        for path in reversed(created):
            path.unlink(missing_ok=True)
        print(f"IQR screening stopped: {exc}", file=sys.stderr)
        print("Verification: FAIL; screening is not complete. Existing outputs were preserved.", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
