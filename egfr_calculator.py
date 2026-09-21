#!/usr/bin/env python3
"""Adult CKD-EPI eGFR calculator with KDIGO G/A categorization."""
from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass, field
from typing import Optional

MGDL_TO_UMOLL = 88.4
AGE_RANGE = (18, 100)
CREATININE_MGDL_RANGE = (0.1, 15.0)
CYSTATIN_C_MGL_RANGE = (0.3, 8.0)
ACR_MGG_RANGE = (0.0, 30000.0)


def _finite(value: float) -> bool:
    return math.isfinite(float(value))


def creatinine_to_mgdl(value: float, unit: str) -> float:
    unit = unit.strip().lower().replace("µ", "u").replace("μ", "u")
    if unit in {"mg/dl", "mgdl", "mg_dl"}:
        return value
    if unit in {"umol/l", "umoll", "umol_l", "micromol/l"}:
        return value / MGDL_TO_UMOLL
    raise ValueError(f"Unsupported creatinine unit: {unit!r} (use mg/dL or µmol/L)")


def _validate(age: float, sex: str, creatinine=None, cystatin=None, acr=None) -> None:
    if not _finite(age) or age < 18:
        raise ValueError("CKD-EPI adult equations require a finite age of at least 18 years.")
    if sex.upper() not in {"M", "F"}:
        raise ValueError("sex must be 'M' or 'F'.")
    if creatinine is not None and (not _finite(creatinine) or creatinine <= 0):
        raise ValueError("Creatinine must be a finite value greater than 0 mg/dL.")
    if cystatin is not None and (not _finite(cystatin) or cystatin <= 0):
        raise ValueError("Cystatin C must be a finite value greater than 0 mg/L.")
    if acr is not None and (not _finite(acr) or acr < 0):
        raise ValueError("ACR must be a finite non-negative value in mg/g.")


def validate_inputs(age: float, sex: str, creatinine_mgdl=None, cystatin_c_mgl=None, acr_mg_g=None) -> list[str]:
    """Return warnings for unusual but computable values; invalid values are rejected by calculations."""
    warnings: list[str] = []
    if sex.upper() not in {"M", "F"}:
        warnings.append(f"Sex must be 'M' or 'F', got {sex!r}.")
    if _finite(age) and age > AGE_RANGE[1]:
        warnings.append(f"Age {age:g} exceeds {AGE_RANGE[1]} years; interpret the estimate with caution.")
    if creatinine_mgdl is not None and _finite(creatinine_mgdl) and not (CREATININE_MGDL_RANGE[0] <= creatinine_mgdl <= CREATININE_MGDL_RANGE[1]):
        warnings.append("Creatinine is outside the calculator's usual input range; verify the value and units.")
    if cystatin_c_mgl is not None and _finite(cystatin_c_mgl) and not (CYSTATIN_C_MGL_RANGE[0] <= cystatin_c_mgl <= CYSTATIN_C_MGL_RANGE[1]):
        warnings.append("Cystatin C is outside the calculator's usual input range; verify the value and units.")
    if acr_mg_g is not None and _finite(acr_mg_g) and acr_mg_g > ACR_MGG_RANGE[1]:
        warnings.append("ACR is extremely high; verify the value and units.")
    return warnings


def ckd_epi_2021_creatinine(age: float, sex: str, scr_mgdl: float) -> float:
    _validate(age, sex, creatinine=scr_mgdl)
    kappa, alpha, factor = (0.7, -0.241, 1.012) if sex.upper() == "F" else (0.9, -0.302, 1.0)
    ratio = scr_mgdl / kappa
    return 142 * min(ratio, 1) ** alpha * max(ratio, 1) ** -1.2 * 0.9938**age * factor


def ckd_epi_2012_cystatin(age: float, sex: str, scys_mgl: float) -> float:
    _validate(age, sex, cystatin=scys_mgl)
    ratio = scys_mgl / 0.8
    factor = 0.932 if sex.upper() == "F" else 1.0
    return 133 * min(ratio, 1) ** -0.499 * max(ratio, 1) ** -1.328 * 0.996**age * factor


def ckd_epi_2021_cystatin(age: float, sex: str, scys_mgl: float) -> float:
    """Compatibility alias for the correctly identified 2012 cystatin-C equation."""
    return ckd_epi_2012_cystatin(age, sex, scys_mgl)


def ckd_epi_2021_creatinine_cystatin(age: float, sex: str, scr_mgdl: float, scys_mgl: float) -> float:
    _validate(age, sex, creatinine=scr_mgdl, cystatin=scys_mgl)
    kappa, alpha, factor = (0.7, -0.219, 0.963) if sex.upper() == "F" else (0.9, -0.144, 1.0)
    cr, cy = scr_mgdl / kappa, scys_mgl / 0.8
    return 135 * min(cr, 1) ** alpha * max(cr, 1) ** -0.544 * min(cy, 1) ** -0.323 * max(cy, 1) ** -0.778 * 0.9961**age * factor


def ckd_epi_2009_creatinine(age: float, sex: str, scr_mgdl: float, race: str = "other") -> float:
    _validate(age, sex, creatinine=scr_mgdl)
    race = race.strip().lower()
    if race not in {"black", "other"}:
        raise ValueError("race must be 'black' or 'other' for the legacy 2009 equation.")
    rf = 1.159 if race == "black" else 1.0
    if sex.upper() == "F":
        exponent = -0.329 if scr_mgdl <= 0.7 else -1.209
        return 144 * (scr_mgdl / 0.7) ** exponent * 0.993**age * rf
    exponent = -0.411 if scr_mgdl <= 0.9 else -1.209
    return 141 * (scr_mgdl / 0.9) ** exponent * 0.993**age * rf


def classify_gfr_stage(egfr: float) -> tuple[str, str]:
    if not _finite(egfr) or egfr < 0:
        raise ValueError("eGFR must be a finite non-negative value.")
    for threshold, label, text in [
        (90, "G1", "Normal or high"), (60, "G2", "Mildly decreased"),
        (45, "G3a", "Mildly to moderately decreased"), (30, "G3b", "Moderately to severely decreased"),
        (15, "G4", "Severely decreased"), (0, "G5", "Kidney failure")]:
        if egfr >= threshold:
            return label, text
    raise AssertionError("unreachable")


def classify_albuminuria(acr_mg_g: float) -> tuple[str, str]:
    if not _finite(acr_mg_g) or acr_mg_g < 0:
        raise ValueError("ACR must be a finite non-negative value.")
    if acr_mg_g < 30:
        return "A1", "Normal to mildly increased"
    if acr_mg_g < 300:
        return "A2", "Moderately increased"
    return "A3", "Severely increased"


KDIGO_RISK_GRID = {
    "G1": {"A1": "Low risk", "A2": "Moderately increased risk", "A3": "High risk"},
    "G2": {"A1": "Low risk", "A2": "Moderately increased risk", "A3": "High risk"},
    "G3a": {"A1": "Moderately increased risk", "A2": "High risk", "A3": "Very high risk"},
    "G3b": {"A1": "High risk", "A2": "Very high risk", "A3": "Very high risk"},
    "G4": {"A1": "Very high risk", "A2": "Very high risk", "A3": "Very high risk"},
    "G5": {"A1": "Very high risk", "A2": "Very high risk", "A3": "Very high risk"},
}


def kdigo_combined_risk(gfr_stage: str, albuminuria_stage: str) -> str:
    try:
        return KDIGO_RISK_GRID[gfr_stage][albuminuria_stage]
    except KeyError as exc:
        raise ValueError("Unknown KDIGO category pair.") from exc


@dataclass
class EgfrResult:
    patient_id: str
    age: float
    sex: str
    egfr_2021_creatinine: Optional[float] = None
    egfr_2012_cystatin: Optional[float] = None
    egfr_2021_combined: Optional[float] = None
    egfr_2009_creatinine: Optional[float] = None
    gfr_stage: Optional[str] = None
    gfr_stage_description: Optional[str] = None
    albuminuria_stage: Optional[str] = None
    albuminuria_stage_description: Optional[str] = None
    combined_risk: Optional[str] = None
    warnings: list[str] = field(default_factory=list)

    @property
    def egfr_2021_cystatin(self) -> Optional[float]:
        return self.egfr_2012_cystatin

    def primary_egfr(self) -> Optional[float]:
        return self.egfr_2021_combined if self.egfr_2021_combined is not None else self.egfr_2012_cystatin if self.egfr_2012_cystatin is not None else self.egfr_2021_creatinine


def calculate_patient(patient_id: str, age: float, sex: str, creatinine_mgdl=None, cystatin_c_mgl=None, acr_mg_g=None, race: str = "other") -> EgfrResult:
    _validate(age, sex, creatinine_mgdl, cystatin_c_mgl, acr_mg_g)
    if creatinine_mgdl is None and cystatin_c_mgl is None:
        raise ValueError("At least one of creatinine or cystatin C is required.")
    r = EgfrResult(patient_id, age, sex.upper(), warnings=validate_inputs(age, sex, creatinine_mgdl, cystatin_c_mgl, acr_mg_g))
    if creatinine_mgdl is not None:
        r.egfr_2021_creatinine = ckd_epi_2021_creatinine(age, sex, creatinine_mgdl)
        r.egfr_2009_creatinine = ckd_epi_2009_creatinine(age, sex, creatinine_mgdl, race)
    if cystatin_c_mgl is not None:
        r.egfr_2012_cystatin = ckd_epi_2012_cystatin(age, sex, cystatin_c_mgl)
    if creatinine_mgdl is not None and cystatin_c_mgl is not None:
        r.egfr_2021_combined = ckd_epi_2021_creatinine_cystatin(age, sex, creatinine_mgdl, cystatin_c_mgl)
    r.gfr_stage, r.gfr_stage_description = classify_gfr_stage(r.primary_egfr())
    if acr_mg_g is not None:
        r.albuminuria_stage, r.albuminuria_stage_description = classify_albuminuria(acr_mg_g)
        r.combined_risk = kdigo_combined_risk(r.gfr_stage, r.albuminuria_stage)
    return r


CSV_INPUT_FIELDS = ["patient_id", "age", "sex", "creatinine", "creatinine_unit", "cystatin_c", "acr", "race"]
CSV_OUTPUT_FIELDS = ["patient_id", "age", "sex", "egfr_2021_creatinine", "egfr_2012_cystatin", "egfr_2021_combined", "egfr_2009_creatinine", "primary_egfr", "gfr_category", "gfr_category_description", "albuminuria_category", "albuminuria_category_description", "kdigo_risk", "warnings"]


def _optional_float(value, label: str, warnings: list[str]):
    if value is None or not str(value).strip():
        return None
    try:
        number = float(value)
        if not math.isfinite(number):
            raise ValueError
        return number
    except (TypeError, ValueError):
        warnings.append(f"Could not parse {label} value {value!r}; ignored.")
        return None


def _fmt(value):
    return "" if value is None else f"{value:.2f}"


def process_csv(input_path: str, output_path: str) -> list[EgfrResult]:
    results: list[EgfrResult] = []
    with open(input_path, newline="", encoding="utf-8-sig") as src:
        reader = csv.DictReader(src)
        missing = {"patient_id", "age", "sex"} - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Input CSV is missing required columns: {sorted(missing)}")
        for n, row in enumerate(reader, 2):
            pid = (row.get("patient_id") or "").strip() or f"row{n}"
            warnings: list[str] = []
            age = _optional_float(row.get("age"), "age", warnings)
            sex = (row.get("sex") or "").strip().upper()
            if age is None or age < 18 or sex not in {"M", "F"}:
                if age is not None and age < 18:
                    warnings.append("CKD-EPI adult equations require age >=18 years.")
                if sex not in {"M", "F"}:
                    warnings.append(f"Sex must be 'M' or 'F', got {row.get('sex')!r}.")
                results.append(EgfrResult(pid, age if age is not None else float("nan"), sex or "?", warnings=warnings))
                continue
            cr_raw = _optional_float(row.get("creatinine"), "creatinine", warnings)
            cystatin = _optional_float(row.get("cystatin_c"), "cystatin C", warnings)
            acr = _optional_float(row.get("acr"), "ACR", warnings)
            cr = None
            if cr_raw is not None:
                try:
                    cr = creatinine_to_mgdl(cr_raw, (row.get("creatinine_unit") or "mg/dL").strip() or "mg/dL")
                    if cr <= 0:
                        raise ValueError("Creatinine must be greater than 0.")
                except ValueError as exc:
                    warnings.append(str(exc)); cr = None
            if cystatin is not None and cystatin <= 0:
                warnings.append("Cystatin C must be greater than 0."); cystatin = None
            if acr is not None and acr < 0:
                warnings.append("ACR must be non-negative."); acr = None
            if cr is None and cystatin is None:
                warnings.append("No valid creatinine or cystatin C provided; eGFR was not computed.")
                results.append(EgfrResult(pid, age, sex, warnings=warnings)); continue
            race = (row.get("race") or "other").strip().lower() or "other"
            if race not in {"black", "other"}:
                warnings.append(f"Unrecognized legacy race value {race!r}; using 'other' for the 2009 comparison."); race = "other"
            try:
                result = calculate_patient(pid, age, sex, cr, cystatin, acr, race)
            except ValueError as exc:
                warnings.append(str(exc)); result = EgfrResult(pid, age, sex)
            result.warnings = warnings + result.warnings
            results.append(result)
    with open(output_path, "w", newline="", encoding="utf-8") as dst:
        writer = csv.DictWriter(dst, fieldnames=CSV_OUTPUT_FIELDS); writer.writeheader()
        for r in results:
            writer.writerow({
                "patient_id": r.patient_id, "age": r.age, "sex": r.sex,
                "egfr_2021_creatinine": _fmt(r.egfr_2021_creatinine), "egfr_2012_cystatin": _fmt(r.egfr_2012_cystatin),
                "egfr_2021_combined": _fmt(r.egfr_2021_combined), "egfr_2009_creatinine": _fmt(r.egfr_2009_creatinine),
                "primary_egfr": _fmt(r.primary_egfr()), "gfr_category": r.gfr_stage or "", "gfr_category_description": r.gfr_stage_description or "",
                "albuminuria_category": r.albuminuria_stage or "", "albuminuria_category_description": r.albuminuria_stage_description or "",
                "kdigo_risk": r.combined_risk or "", "warnings": " | ".join(r.warnings)})
    return results


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="egfr_calculator", description="Adult CKD-EPI eGFR calculator with KDIGO G/A categorization.")
    sub = parser.add_subparsers(dest="command", required=True)
    single = sub.add_parser("single"); single.add_argument("--id", dest="patient_id", default="patient"); single.add_argument("--age", type=float, required=True)
    single.add_argument("--sex", choices=["M", "F", "m", "f"], required=True); single.add_argument("--creatinine", type=float); single.add_argument("--creatinine-unit", default="mg/dL")
    single.add_argument("--cystatin-c", type=float); single.add_argument("--acr", type=float); single.add_argument("--race", default="other", choices=["black", "other"])
    batch = sub.add_parser("batch"); batch.add_argument("-i", "--input", required=True); batch.add_argument("-o", "--output", required=True)
    return parser


def _print_single_result(r: EgfrResult) -> None:
    print(f"Patient: {r.patient_id}  Age: {r.age:g}  Sex: {r.sex}")
    for label, value in [("eGFRcr (2021 CKD-EPI)", r.egfr_2021_creatinine), ("eGFRcys (2012 CKD-EPI)", r.egfr_2012_cystatin), ("eGFRcr-cys (2021 CKD-EPI)", r.egfr_2021_combined), ("eGFRcr (2009 CKD-EPI, legacy)", r.egfr_2009_creatinine)]:
        if value is not None: print(f"  {label}: {value:.1f} mL/min/1.73 m²")
    if r.primary_egfr() is not None:
        print(f"  Primary eGFR: {r.primary_egfr():.1f} mL/min/1.73 m²")
        print(f"  GFR category: {r.gfr_stage} ({r.gfr_stage_description})")
        if r.albuminuria_stage:
            print(f"  Albuminuria category: {r.albuminuria_stage} ({r.albuminuria_stage_description})")
            print(f"  KDIGO risk category: {r.combined_risk}")
    for warning in r.warnings: print(f"  Warning: {warning}")
    print("  Note: a single result does not establish chronic kidney disease.")


def main(argv: Optional[list[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    if args.command == "single":
        if args.creatinine is None and args.cystatin_c is None:
            print("Error: at least one of --creatinine or --cystatin-c is required.", file=sys.stderr); return 2
        try:
            cr = None if args.creatinine is None else creatinine_to_mgdl(args.creatinine, args.creatinine_unit)
            result = calculate_patient(args.patient_id, args.age, args.sex, cr, args.cystatin_c, args.acr, args.race)
        except ValueError as exc:
            print(f"Error: {exc}", file=sys.stderr); return 2
        _print_single_result(result); return 0
    try:
        results = process_csv(args.input, args.output)
    except (OSError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr); return 2
    print(f"Processed {len(results)} rows -> {args.output}")
    warned = sum(bool(r.warnings) for r in results)
    if warned: print(f"{warned} row(s) contain warnings; review the output CSV.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
