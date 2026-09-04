#!/usr/bin/env python3
"""
CKD-EPI eGFR Calculator
========================

Computes estimated glomerular filtration rate (eGFR) using:
  - 2021 CKD-EPI creatinine equation (race-free, Inker et al., NEJM 2021)
  - 2021 CKD-EPI cystatin C equation (race-free)
  - 2021 CKD-EPI creatinine-cystatin C combined equation (race-free)
  - 2009 CKD-EPI creatinine equation (with race coefficient, Levey et al.)

Stages chronic kidney disease per KDIGO 2012 GFR/albuminuria grid.

Stdlib only. Usage: python egfr_calculator.py --help
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass, field
from typing import Optional


# ---------------------------------------------------------------------------
# Unit conversion
# ---------------------------------------------------------------------------

MGDL_TO_UMOLL = 88.4  # 1 mg/dL creatinine = 88.4 umol/L


def creatinine_to_mgdl(value: float, unit: str) -> float:
    """Convert a serum creatinine value to mg/dL."""
    unit = unit.strip().lower()
    if unit in ("mg/dl", "mgdl", "mg_dl"):
        return value
    if unit in ("umol/l", "umoll", "umol_l", "micromol/l"):
        return value / MGDL_TO_UMOLL
    raise ValueError(f"Unsupported creatinine unit: {unit!r} (use mg/dL or umol/L)")


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

AGE_RANGE = (18, 100)
CREATININE_MGDL_RANGE = (0.1, 15.0)
CYSTATIN_C_MGL_RANGE = (0.3, 8.0)
ACR_MGG_RANGE = (0.0, 30000.0)


def validate_inputs(
    age: float,
    sex: str,
    creatinine_mgdl: Optional[float] = None,
    cystatin_c_mgl: Optional[float] = None,
    acr_mg_g: Optional[float] = None,
) -> list[str]:
    """Return a list of human-readable warnings for implausible/out-of-range inputs."""
    warnings: list[str] = []

    if sex.upper() not in ("M", "F"):
        warnings.append(f"Sex must be 'M' or 'F', got {sex!r}.")

    if age is None or age != age:
        warnings.append("Age is missing or not a number.")
    else:
        if age < 0:
            warnings.append(f"Age {age} is negative and physiologically impossible.")
        elif age < AGE_RANGE[0]:
            warnings.append(
                f"Age {age} is below {AGE_RANGE[0]}: CKD-EPI equations are "
                f"validated for adults; pediatric equations should be used instead."
            )
        elif age > AGE_RANGE[1]:
            warnings.append(
                f"Age {age} exceeds {AGE_RANGE[1]}: outside the typical validated "
                f"range, interpret result with caution."
            )

    if creatinine_mgdl is not None:
        if creatinine_mgdl <= 0:
            warnings.append(
                f"Creatinine {creatinine_mgdl:.3f} mg/dL is non-positive and implausible."
            )
        elif not (CREATININE_MGDL_RANGE[0] <= creatinine_mgdl <= CREATININE_MGDL_RANGE[1]):
            warnings.append(
                f"Creatinine {creatinine_mgdl:.3f} mg/dL is outside the plausible "
                f"physiologic range {CREATININE_MGDL_RANGE[0]}-{CREATININE_MGDL_RANGE[1]} mg/dL."
            )

    if cystatin_c_mgl is not None:
        if cystatin_c_mgl <= 0:
            warnings.append(
                f"Cystatin C {cystatin_c_mgl:.3f} mg/L is non-positive and implausible."
            )
        elif not (CYSTATIN_C_MGL_RANGE[0] <= cystatin_c_mgl <= CYSTATIN_C_MGL_RANGE[1]):
            warnings.append(
                f"Cystatin C {cystatin_c_mgl:.3f} mg/L is outside the plausible "
                f"physiologic range {CYSTATIN_C_MGL_RANGE[0]}-{CYSTATIN_C_MGL_RANGE[1]} mg/L."
            )

    if acr_mg_g is not None:
        if acr_mg_g < 0:
            warnings.append(f"ACR {acr_mg_g:.1f} mg/g is negative and implausible.")
        elif acr_mg_g > ACR_MGG_RANGE[1]:
            warnings.append(
                f"ACR {acr_mg_g:.1f} mg/g exceeds {ACR_MGG_RANGE[1]} mg/g, an "
                f"extreme value: verify the measurement."
            )

    return warnings


# ---------------------------------------------------------------------------
# 2021 CKD-EPI equations (Inker et al., NEJM 2021) -- race-free
# ---------------------------------------------------------------------------


def ckd_epi_2021_creatinine(age: float, sex: str, scr_mgdl: float) -> float:
    """2021 CKD-EPI creatinine-only equation (race-free).

    eGFR = 142 × min(Scr/κ, 1)^α × max(Scr/κ, 1)^(-1.200) × 0.9938^Age × sex_factor

    Where:
        κ = 0.7 (female) or 0.9 (male)
        α = -0.241 (female) or -0.302 (male)
        sex_factor = 1.012 (female) or 1.0 (male)

    Returns eGFR in mL/min/1.73m².
    """
    sex = sex.upper()
    if sex == "F":
        kappa, alpha, sex_factor = 0.7, -0.241, 1.012
    elif sex == "M":
        kappa, alpha, sex_factor = 0.9, -0.302, 1.0
    else:
        raise ValueError(f"sex must be 'M' or 'F', got {sex!r}")

    ratio = scr_mgdl / kappa
    egfr = (
        142
        * min(ratio, 1.0) ** alpha
        * max(ratio, 1.0) ** -1.200
        * (0.9938 ** age)
        * sex_factor
    )
    return egfr


def ckd_epi_2009_creatinine(age: float, sex: str, scr_mgdl: float, race: str = "other") -> float:
    """2009 CKD-EPI creatinine equation (Levey et al., Ann Intern Med 2009).

    Includes the race coefficient that was removed in the 2021 equation.

    For females, Scr ≤ 0.7:
        eGFR = 144 × (Scr/0.7)^(-0.329) × 0.993^Age × (1.159 if Black)
    For females, Scr > 0.7:
        eGFR = 144 × (Scr/0.7)^(-1.209) × 0.993^Age × (1.159 if Black)
    For males, Scr ≤ 0.9:
        eGFR = 141 × (Scr/0.9)^(-0.411) × 0.993^Age × (1.159 if Black)
    For males, Scr > 0.9:
        eGFR = 141 × (Scr/0.9)^(-1.209) × 0.993^Age × (1.159 if Black)

    Returns eGFR in mL/min/1.73m².
    """
    sex = sex.upper()
    race = race.lower()
    race_factor = 1.159 if race == "black" else 1.0

    if sex == "F":
        if scr_mgdl <= 0.7:
            egfr = 144 * (scr_mgdl / 0.7) ** (-0.329) * 0.993 ** age * race_factor
        else:
            egfr = 144 * (scr_mgdl / 0.7) ** (-1.209) * 0.993 ** age * race_factor
    elif sex == "M":
        if scr_mgdl <= 0.9:
            egfr = 141 * (scr_mgdl / 0.9) ** (-0.411) * 0.993 ** age * race_factor
        else:
            egfr = 141 * (scr_mgdl / 0.9) ** (-1.209) * 0.993 ** age * race_factor
    else:
        raise ValueError(f"sex must be 'M' or 'F', got {sex!r}")

    return egfr


def ckd_epi_2021_cystatin(age: float, sex: str, scys_mgl: float) -> float:
    """2021 CKD-EPI cystatin C-only equation (race-free).

    Returns eGFR in mL/min/1.73m².
    """
    sex = sex.upper()
    if sex == "F":
        sex_factor = 0.932
    elif sex == "M":
        sex_factor = 1.0
    else:
        raise ValueError(f"sex must be 'M' or 'F', got {sex!r}")

    ratio = scys_mgl / 0.8
    egfr = (
        133
        * min(ratio, 1.0) ** -0.499
        * max(ratio, 1.0) ** -1.328
        * (0.9961 ** age)
        * sex_factor
    )
    return egfr


def ckd_epi_2021_creatinine_cystatin(
    age: float, sex: str, scr_mgdl: float, scys_mgl: float
) -> float:
    """2021 CKD-EPI creatinine-cystatin C combined equation (race-free).

    Returns eGFR in mL/min/1.73m².
    """
    sex = sex.upper()
    if sex == "F":
        kappa, alpha, sex_factor = 0.7, -0.219, 0.963
    elif sex == "M":
        kappa, alpha, sex_factor = 0.9, -0.144, 1.0
    else:
        raise ValueError(f"sex must be 'M' or 'F', got {sex!r}")

    scr_ratio = scr_mgdl / kappa
    scys_ratio = scys_mgl / 0.8
    egfr = (
        135
        * min(scr_ratio, 1.0) ** alpha
        * max(scr_ratio, 1.0) ** -0.544
        * min(scys_ratio, 1.0) ** -0.323
        * max(scys_ratio, 1.0) ** -0.778
        * (0.9961 ** age)
        * sex_factor
    )
    return egfr


# ---------------------------------------------------------------------------
# KDIGO 2012 CKD staging
# ---------------------------------------------------------------------------

GFR_STAGES = (
    ("G1", 90.0, float("inf"), "Normal or high"),
    ("G2", 60.0, 89.999999, "Mildly decreased"),
    ("G3a", 45.0, 59.999999, "Mildly to moderately decreased"),
    ("G3b", 30.0, 44.999999, "Moderately to severely decreased"),
    ("G4", 15.0, 29.999999, "Severely decreased"),
    ("G5", 0.0, 14.999999, "Kidney failure"),
)

ALBUMINURIA_STAGES = (
    ("A1", 0.0, 29.999999, "Normal to mildly increased"),
    ("A2", 30.0, 299.999999, "Moderately increased"),
    ("A3", 300.0, float("inf"), "Severely increased"),
)

KDIGO_RISK_GRID = {
    ("G1", "A1"): "Low risk",
    ("G1", "A2"): "Moderately increased risk",
    ("G1", "A3"): "High risk",
    ("G2", "A1"): "Low risk",
    ("G2", "A2"): "Moderately increased risk",
    ("G2", "A3"): "High risk",
    ("G3a", "A1"): "Moderately increased risk",
    ("G3a", "A2"): "High risk",
    ("G3a", "A3"): "Very high risk",
    ("G3b", "A1"): "High risk",
    ("G3b", "A2"): "Very high risk",
    ("G3b", "A3"): "Very high risk",
    ("G4", "A1"): "Very high risk",
    ("G4", "A2"): "Very high risk",
    ("G4", "A3"): "Very high risk",
    ("G5", "A1"): "Very high risk",
    ("G5", "A2"): "Very high risk",
    ("G5", "A3"): "Very high risk",
}


def classify_gfr_stage(egfr: float) -> tuple[str, str]:
    """Return (stage_label, description) for a KDIGO GFR category."""
    for label, low, high, description in GFR_STAGES:
        if low <= egfr <= high:
            return label, description
    return "G5", "Kidney failure"


def classify_albuminuria(acr_mg_g: float) -> tuple[str, str]:
    """Return (stage_label, description) for a KDIGO albuminuria category."""
    for label, low, high, description in ALBUMINURIA_STAGES:
        if low <= acr_mg_g <= high:
            return label, description
    return "A3", "Severely increased"


def kdigo_combined_risk(gfr_stage: str, albuminuria_stage: str) -> str:
    """Look up the KDIGO 2012 heat-map risk category for a G/A stage pair."""
    return KDIGO_RISK_GRID.get((gfr_stage, albuminuria_stage), "Unknown")


# ---------------------------------------------------------------------------
# Patient-level result assembly
# ---------------------------------------------------------------------------


@dataclass
class EgfrResult:
    patient_id: str
    age: float
    sex: str
    egfr_2021_creatinine: Optional[float] = None
    egfr_2009_creatinine: Optional[float] = None
    egfr_2021_cystatin: Optional[float] = None
    egfr_2021_combined: Optional[float] = None
    gfr_stage: Optional[str] = None
    gfr_stage_description: Optional[str] = None
    albuminuria_stage: Optional[str] = None
    albuminuria_stage_description: Optional[str] = None
    combined_risk: Optional[str] = None
    warnings: list[str] = field(default_factory=list)

    def primary_egfr(self) -> Optional[float]:
        """Best available eGFR, preferring combined > cystatin > 2021 creatinine."""
        if self.egfr_2021_combined is not None:
            return self.egfr_2021_combined
        if self.egfr_2021_cystatin is not None:
            return self.egfr_2021_cystatin
        return self.egfr_2021_creatinine


def calculate_patient(
    patient_id: str,
    age: float,
    sex: str,
    creatinine_mgdl: Optional[float] = None,
    cystatin_c_mgl: Optional[float] = None,
    acr_mg_g: Optional[float] = None,
    race: str = "other",
) -> EgfrResult:
    """Run all applicable equations and KDIGO staging for one patient."""
    warnings = validate_inputs(age, sex, creatinine_mgdl, cystatin_c_mgl, acr_mg_g)
    result = EgfrResult(patient_id=patient_id, age=age, sex=sex.upper(), warnings=warnings)

    if creatinine_mgdl is not None:
        result.egfr_2021_creatinine = ckd_epi_2021_creatinine(age, sex, creatinine_mgdl)
        result.egfr_2009_creatinine = ckd_epi_2009_creatinine(age, sex, creatinine_mgdl, race)

    if cystatin_c_mgl is not None:
        result.egfr_2021_cystatin = ckd_epi_2021_cystatin(age, sex, cystatin_c_mgl)

    if creatinine_mgdl is not None and cystatin_c_mgl is not None:
        result.egfr_2021_combined = ckd_epi_2021_creatinine_cystatin(
            age, sex, creatinine_mgdl, cystatin_c_mgl
        )

    primary = result.primary_egfr()
    if primary is not None:
        stage, description = classify_gfr_stage(primary)
        result.gfr_stage = stage
        result.gfr_stage_description = description

        if acr_mg_g is not None:
            a_stage, a_description = classify_albuminuria(acr_mg_g)
            result.albuminuria_stage = a_stage
            result.albuminuria_stage_description = a_description
            result.combined_risk = kdigo_combined_risk(stage, a_stage)

    return result


# ---------------------------------------------------------------------------
# CSV batch processing
# ---------------------------------------------------------------------------

CSV_INPUT_FIELDS = [
    "patient_id", "age", "sex", "creatinine", "creatinine_unit",
    "cystatin_c", "acr", "race",
]

CSV_OUTPUT_FIELDS = [
    "patient_id", "age", "sex",
    "egfr_2021_creatinine", "egfr_2009_creatinine",
    "egfr_2021_cystatin", "egfr_2021_combined",
    "primary_egfr", "gfr_stage", "gfr_stage_description",
    "albuminuria_stage", "albuminuria_stage_description",
    "combined_risk", "warnings",
]


def _parse_optional_float(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    value = value.strip()
    if value == "":
        return None
    return float(value)


def process_csv(input_path: str, output_path: str) -> list[EgfrResult]:
    """Read patient rows from a CSV, compute eGFR/staging for each, write results CSV."""
    results: list[EgfrResult] = []

    with open(input_path, "r", newline="", encoding="utf-8-sig") as f_in:
        reader = csv.DictReader(f_in)
        missing = set(["patient_id", "age", "sex"]) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Input CSV is missing required columns: {sorted(missing)}")

        for row_num, row in enumerate(reader, start=2):
            patient_id = (row.get("patient_id") or "").strip() or f"row{row_num}"
            row_warnings: list[str] = []

            try:
                age = float(row["age"])
            except (KeyError, ValueError, TypeError):
                row_warnings.append(f"Could not parse age {row.get('age')!r}; row skipped.")
                results.append(EgfrResult(patient_id=patient_id, age=float("nan"), sex="?",
                                           warnings=row_warnings))
                continue

            sex = (row.get("sex") or "").strip().upper()
            race = (row.get("race") or "other").strip().lower()

            creatinine_raw = _parse_optional_float(row.get("creatinine"))
            creatinine_unit = (row.get("creatinine_unit") or "mg/dL").strip() or "mg/dL"
            creatinine_mgdl = None
            if creatinine_raw is not None:
                try:
                    creatinine_mgdl = creatinine_to_mgdl(creatinine_raw, creatinine_unit)
                except ValueError as exc:
                    row_warnings.append(str(exc))

            cystatin_c_mgl = _parse_optional_float(row.get("cystatin_c"))
            acr_mg_g = _parse_optional_float(row.get("acr"))

            if sex not in ("M", "F"):
                row_warnings.append(f"Sex must be 'M' or 'F', got {row.get('sex')!r}.")
                results.append(EgfrResult(patient_id=patient_id, age=age, sex=sex or "?",
                                           warnings=row_warnings))
                continue

            if creatinine_mgdl is None and cystatin_c_mgl is None:
                row_warnings.append("No creatinine or cystatin C provided; cannot compute eGFR.")
                results.append(EgfrResult(patient_id=patient_id, age=age, sex=sex,
                                           warnings=row_warnings))
                continue

            result = calculate_patient(
                patient_id=patient_id, age=age, sex=sex,
                creatinine_mgdl=creatinine_mgdl, cystatin_c_mgl=cystatin_c_mgl,
                acr_mg_g=acr_mg_g, race=race,
            )
            result.warnings = row_warnings + result.warnings
            results.append(result)

    with open(output_path, "w", newline="", encoding="utf-8") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=CSV_OUTPUT_FIELDS)
        writer.writeheader()
        for r in results:
            writer.writerow({
                "patient_id": r.patient_id,
                "age": r.age,
                "sex": r.sex,
                "egfr_2021_creatinine": _fmt(r.egfr_2021_creatinine),
                "egfr_2009_creatinine": _fmt(r.egfr_2009_creatinine),
                "egfr_2021_cystatin": _fmt(r.egfr_2021_cystatin),
                "egfr_2021_combined": _fmt(r.egfr_2021_combined),
                "primary_egfr": _fmt(r.primary_egfr()),
                "gfr_stage": r.gfr_stage or "",
                "gfr_stage_description": r.gfr_stage_description or "",
                "albuminuria_stage": r.albuminuria_stage or "",
                "albuminuria_stage_description": r.albuminuria_stage_description or "",
                "combined_risk": r.combined_risk or "",
                "warnings": " | ".join(r.warnings),
            })

    return results


def _fmt(value: Optional[float]) -> str:
    return "" if value is None else f"{value:.2f}"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="egfr_calculator",
        description=(
            "CKD-EPI eGFR calculator with KDIGO 2012 CKD staging. "
            "Computes eGFR from creatinine and/or cystatin C."
        ),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    single = subparsers.add_parser("single", help="Calculate eGFR for one patient")
    single.add_argument("--id", dest="patient_id", default="patient", help="Patient identifier")
    single.add_argument("--age", type=float, required=True, help="Age in years")
    single.add_argument("--sex", required=True, choices=["M", "F", "m", "f"], help="Biological sex (M/F)")
    single.add_argument("--creatinine", type=float, default=None, help="Serum creatinine value")
    single.add_argument("--creatinine-unit", default="mg/dL",
                        choices=["mg/dL", "mg/dl", "umol/L", "umol/l"],
                        help="Unit of the creatinine value (default: mg/dL)")
    single.add_argument("--cystatin-c", type=float, default=None, help="Serum cystatin C in mg/L")
    single.add_argument("--acr", type=float, default=None, help="Urine albumin-creatinine ratio in mg/g")
    single.add_argument("--race", default="other", choices=["black", "other"],
                        help="Race for 2009 CKD-EPI equation (default: other)")

    batch = subparsers.add_parser("batch", help="Calculate eGFR for a CSV of patients")
    batch.add_argument("-i", "--input", required=True, help="Path to input CSV")
    batch.add_argument("-o", "--output", required=True, help="Path to write output CSV")

    return parser


def _print_single_result(result: EgfrResult) -> None:
    print(f"Patient: {result.patient_id}  Age: {result.age}  Sex: {result.sex}")
    if result.egfr_2021_creatinine is not None:
        print(f"  eGFR (2021 CKD-EPI creatinine, race-free):    {result.egfr_2021_creatinine:.1f} mL/min/1.73m²")
    if result.egfr_2009_creatinine is not None:
        print(f"  eGFR (2009 CKD-EPI creatinine):               {result.egfr_2009_creatinine:.1f} mL/min/1.73m²")
    if result.egfr_2021_cystatin is not None:
        print(f"  eGFR (2021 CKD-EPI cystatin C, race-free):    {result.egfr_2021_cystatin:.1f} mL/min/1.73m²")
    if result.egfr_2021_combined is not None:
        print(f"  eGFR (2021 CKD-EPI Cr-CysC combined):         {result.egfr_2021_combined:.1f} mL/min/1.73m²")

    primary = result.primary_egfr()
    if primary is None:
        print("  No creatinine or cystatin C provided: eGFR could not be computed.")
    else:
        print(f"  Primary eGFR used for staging: {primary:.1f} mL/min/1.73m²")
        print(f"  KDIGO GFR stage: {result.gfr_stage} ({result.gfr_stage_description})")
        if result.albuminuria_stage:
            print(f"  KDIGO albuminuria stage: {result.albuminuria_stage} ({result.albuminuria_stage_description})")
            print(f"  KDIGO combined risk category: {result.combined_risk}")

    if result.warnings:
        print("  Warnings:")
        for w in result.warnings:
            print(f"    - {w}")


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    if args.command == "single":
        creatinine_mgdl = None
        if args.creatinine is not None:
            try:
                creatinine_mgdl = creatinine_to_mgdl(args.creatinine, args.creatinine_unit)
            except ValueError as exc:
                print(f"Error: {exc}", file=sys.stderr)
                return 1

        if creatinine_mgdl is None and args.cystatin_c is None:
            print("Error: at least one of --creatinine or --cystatin-c is required.", file=sys.stderr)
            return 1

        result = calculate_patient(
            patient_id=args.patient_id, age=args.age, sex=args.sex,
            creatinine_mgdl=creatinine_mgdl, cystatin_c_mgl=args.cystatin_c,
            acr_mg_g=args.acr, race=args.race,
        )
        _print_single_result(result)
        return 0

    if args.command == "batch":
        results = process_csv(args.input, args.output)
        n_warned = sum(1 for r in results if r.warnings)
        print(f"Processed {len(results)} patients -> {args.output}")
        if n_warned:
            print(f"{n_warned} patient(s) had warnings; see the 'warnings' column.")
        return 0

    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
