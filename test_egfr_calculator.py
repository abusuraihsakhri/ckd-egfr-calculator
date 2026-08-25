"""Tests for egfr_calculator.py -- plain assert statements, stdlib only.

Run with: python test_egfr_calculator.py
"""

import csv
import math
import os
import tempfile

import egfr_calculator as egfr


# ---------------------------------------------------------------------------
# CKD-EPI 2021 Creatinine (race-free)
# ---------------------------------------------------------------------------

def test_2021_creatinine_female_above_kappa():
    """50-year-old female, Scr=1.0 (>0.7 kappa), uses -1.200 branch."""
    result = egfr.ckd_epi_2021_creatinine(age=50, sex="F", scr_mgdl=1.0)
    assert math.isclose(result, 68.6, abs_tol=0.5), result


def test_2021_creatinine_male_above_kappa():
    """60-year-old male, Scr=1.2 (>0.9 kappa), uses -1.200 branch."""
    result = egfr.ckd_epi_2021_creatinine(age=60, sex="M", scr_mgdl=1.2)
    assert math.isclose(result, 69.2, abs_tol=0.5), result


def test_2021_creatinine_below_kappa_uses_alpha_branch():
    """Low Scr (below kappa) uses alpha exponent, producing higher eGFR."""
    result = egfr.ckd_epi_2021_creatinine(age=30, sex="F", scr_mgdl=0.5)
    assert result > 90, result


def test_2021_creatinine_male_below_kappa():
    """Male with Scr below 0.9 kappa threshold."""
    result = egfr.ckd_epi_2021_creatinine(age=40, sex="M", scr_mgdl=0.7)
    assert result > 90, result


def test_2021_creatinine_at_kappa_boundary():
    """When Scr == kappa, min and max both equal 1, isolating age/sex terms."""
    # Female: 142 * 1^(-0.241) * 1^(-1.200) * 0.9938^40 * 1.012
    expected = 142 * 1.0 * 1.0 * (0.9938 ** 40) * 1.012
    result = egfr.ckd_epi_2021_creatinine(age=40, sex="F", scr_mgdl=0.7)
    assert math.isclose(result, expected, abs_tol=0.01), (result, expected)


# ---------------------------------------------------------------------------
# CKD-EPI 2009 Creatinine (with race coefficient)
# ---------------------------------------------------------------------------

def test_2009_creatinine_female_nonblack():
    """2009 CKD-EPI: 50yo female, Scr=1.0, non-Black."""
    result = egfr.ckd_epi_2009_creatinine(age=50, sex="F", scr_mgdl=1.0, race="other")
    # Scr > 0.7: 144 * (1.0/0.7)^(-1.209) * 0.993^50 * 1.0
    expected = 144 * (1.0 / 0.7) ** (-1.209) * 0.993 ** 50
    assert math.isclose(result, expected, abs_tol=0.1), (result, expected)


def test_2009_creatinine_male_black():
    """2009 CKD-EPI: 45yo male, Scr=1.0, Black -> race factor 1.159."""
    result = egfr.ckd_epi_2009_creatinine(age=45, sex="M", scr_mgdl=1.0, race="black")
    # Scr > 0.9: 141 * (1.0/0.9)^(-1.209) * 0.993^45 * 1.159
    expected = 141 * (1.0 / 0.9) ** (-1.209) * 0.993 ** 45 * 1.159
    assert math.isclose(result, expected, abs_tol=0.1), (result, expected)


def test_2009_creatinine_below_threshold():
    """2009 CKD-EPI: male Scr=0.7 (<0.9 threshold), uses -0.411 branch."""
    result = egfr.ckd_epi_2009_creatinine(age=30, sex="M", scr_mgdl=0.7, race="other")
    expected = 141 * (0.7 / 0.9) ** (-0.411) * 0.993 ** 30
    assert math.isclose(result, expected, abs_tol=0.1), (result, expected)


def test_2009_vs_2021_nonblack_similar():
    """For non-Black patients, 2009 and 2021 should give similar results."""
    r2009 = egfr.ckd_epi_2009_creatinine(age=50, sex="M", scr_mgdl=1.0, race="other")
    r2021 = egfr.ckd_epi_2021_creatinine(age=50, sex="M", scr_mgdl=1.0)
    # Should be within ~5 mL/min of each other
    assert abs(r2009 - r2021) < 5, (r2009, r2021)


def test_2009_black_higher_than_nonblack():
    """Black race factor should produce higher eGFR in 2009 equation."""
    r_black = egfr.ckd_epi_2009_creatinine(age=50, sex="M", scr_mgdl=1.2, race="black")
    r_other = egfr.ckd_epi_2009_creatinine(age=50, sex="M", scr_mgdl=1.2, race="other")
    assert r_black > r_other, (r_black, r_other)


# ---------------------------------------------------------------------------
# CKD-EPI 2021 Cystatin C
# ---------------------------------------------------------------------------

def test_cystatin_equation_at_threshold():
    """Cystatin C at kappa (0.8) collapses min/max to 1."""
    result = egfr.ckd_epi_2021_cystatin(age=50, sex="F", scys_mgl=0.8)
    expected = 133 * 0.9961 ** 50 * 0.932
    assert math.isclose(result, expected, abs_tol=0.5), result


def test_cystatin_equation_male():
    """Male cystatin C with value above threshold."""
    result = egfr.ckd_epi_2021_cystatin(age=60, sex="M", scys_mgl=1.2)
    expected = 133 * (1.2 / 0.8) ** (-1.328) * 0.9961 ** 60
    assert math.isclose(result, expected, abs_tol=0.5), (result, expected)


# ---------------------------------------------------------------------------
# CKD-EPI 2021 Combined (Cr + CysC)
# ---------------------------------------------------------------------------

def test_combined_equation_between_creatinine_and_cystatin():
    """Combined equation should land between the two single-marker estimates."""
    age, sex, scr, scys = 55, "M", 1.1, 1.0
    cr = egfr.ckd_epi_2021_creatinine(age, sex, scr)
    cy = egfr.ckd_epi_2021_cystatin(age, sex, scys)
    combined = egfr.ckd_epi_2021_creatinine_cystatin(age, sex, scr, scys)
    lo, hi = min(cr, cy), max(cr, cy)
    assert lo - 5 <= combined <= hi + 5, (cr, cy, combined)


# ---------------------------------------------------------------------------
# Sex validation
# ---------------------------------------------------------------------------

def test_sex_must_be_m_or_f():
    try:
        egfr.ckd_epi_2021_creatinine(age=40, sex="X", scr_mgdl=1.0)
        assert False, "expected ValueError for invalid sex"
    except ValueError:
        pass


def test_2009_sex_must_be_m_or_f():
    try:
        egfr.ckd_epi_2009_creatinine(age=40, sex="X", scr_mgdl=1.0)
        assert False, "expected ValueError"
    except ValueError:
        pass


# ---------------------------------------------------------------------------
# GFR Staging
# ---------------------------------------------------------------------------

def test_gfr_staging_g1():
    assert egfr.classify_gfr_stage(120)[0] == "G1"
    assert egfr.classify_gfr_stage(90)[0] == "G1"


def test_gfr_staging_g2():
    assert egfr.classify_gfr_stage(75)[0] == "G2"
    assert egfr.classify_gfr_stage(60)[0] == "G2"


def test_gfr_staging_g3a():
    assert egfr.classify_gfr_stage(50)[0] == "G3a"
    assert egfr.classify_gfr_stage(45)[0] == "G3a"


def test_gfr_staging_g3b():
    assert egfr.classify_gfr_stage(35)[0] == "G3b"
    assert egfr.classify_gfr_stage(30)[0] == "G3b"


def test_gfr_staging_g4():
    assert egfr.classify_gfr_stage(20)[0] == "G4"
    assert egfr.classify_gfr_stage(15)[0] == "G4"


def test_gfr_staging_g5():
    assert egfr.classify_gfr_stage(5)[0] == "G5"
    assert egfr.classify_gfr_stage(0)[0] == "G5"


# ---------------------------------------------------------------------------
# Albuminuria Staging
# ---------------------------------------------------------------------------

def test_albuminuria_a1():
    assert egfr.classify_albuminuria(10)[0] == "A1"
    assert egfr.classify_albuminuria(29.9)[0] == "A1"


def test_albuminuria_a2():
    assert egfr.classify_albuminuria(30)[0] == "A2"
    assert egfr.classify_albuminuria(200)[0] == "A2"


def test_albuminuria_a3():
    assert egfr.classify_albuminuria(300)[0] == "A3"
    assert egfr.classify_albuminuria(1000)[0] == "A3"


# ---------------------------------------------------------------------------
# KDIGO Risk Grid
# ---------------------------------------------------------------------------

def test_kdigo_risk_grid():
    assert egfr.kdigo_combined_risk("G1", "A1") == "Low risk"
    assert egfr.kdigo_combined_risk("G1", "A3") == "High risk"
    assert egfr.kdigo_combined_risk("G3a", "A1") == "Moderately increased risk"
    assert egfr.kdigo_combined_risk("G4", "A1") == "Very high risk"
    assert egfr.kdigo_combined_risk("G5", "A3") == "Very high risk"


# ---------------------------------------------------------------------------
# Unit conversion
# ---------------------------------------------------------------------------

def test_creatinine_unit_conversion_roundtrip():
    mgdl = 1.5
    umoll = mgdl * egfr.MGDL_TO_UMOLL
    assert math.isclose(egfr.creatinine_to_mgdl(umoll, "umol/L"), mgdl, abs_tol=1e-9)
    assert math.isclose(egfr.creatinine_to_mgdl(mgdl, "mg/dL"), mgdl, abs_tol=1e-9)


def test_unsupported_unit_raises():
    try:
        egfr.creatinine_to_mgdl(1.0, "mg/L")
        assert False, "expected ValueError for unsupported unit"
    except ValueError:
        pass


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def test_validate_inputs_flags_implausible_values():
    warnings = egfr.validate_inputs(
        age=-5, sex="F", creatinine_mgdl=50.0, cystatin_c_mgl=0.01
    )
    assert any("negative" in w for w in warnings)
    assert any("Creatinine" in w for w in warnings)
    assert any("Cystatin C" in w for w in warnings)


def test_validate_inputs_accepts_normal_values():
    warnings = egfr.validate_inputs(age=45, sex="M", creatinine_mgdl=1.0, cystatin_c_mgl=0.9)
    assert warnings == [], warnings


# ---------------------------------------------------------------------------
# Patient workflow
# ---------------------------------------------------------------------------

def test_calculate_patient_full_workflow():
    result = egfr.calculate_patient(
        patient_id="P001", age=65, sex="F",
        creatinine_mgdl=1.3, cystatin_c_mgl=1.1, acr_mg_g=150,
    )
    assert result.egfr_2021_creatinine is not None
    assert result.egfr_2009_creatinine is not None
    assert result.egfr_2021_cystatin is not None
    assert result.egfr_2021_combined is not None
    assert result.gfr_stage is not None
    assert result.albuminuria_stage == "A2"
    assert result.combined_risk is not None


def test_calculate_patient_creatinine_only():
    result = egfr.calculate_patient(patient_id="P002", age=40, sex="M", creatinine_mgdl=0.9)
    assert result.egfr_2021_creatinine is not None
    assert result.egfr_2009_creatinine is not None
    assert result.egfr_2021_cystatin is None
    assert result.egfr_2021_combined is None
    assert result.primary_egfr() == result.egfr_2021_creatinine


def test_calculate_patient_with_race():
    """2009 equation should differ by race; 2021 should be identical."""
    r_black = egfr.calculate_patient("PB", 50, "M", creatinine_mgdl=1.2, race="black")
    r_other = egfr.calculate_patient("PO", 50, "M", creatinine_mgdl=1.2, race="other")
    assert r_black.egfr_2021_creatinine == r_other.egfr_2021_creatinine
    assert r_black.egfr_2009_creatinine != r_other.egfr_2009_creatinine


# ---------------------------------------------------------------------------
# CSV batch processing
# ---------------------------------------------------------------------------

def test_batch_csv_processing():
    with tempfile.TemporaryDirectory() as tmp:
        input_path = os.path.join(tmp, "patients.csv")
        output_path = os.path.join(tmp, "results.csv")

        rows = [
            {"patient_id": "A1", "age": "50", "sex": "F", "creatinine": "1.0",
             "creatinine_unit": "mg/dL", "cystatin_c": "", "acr": "", "race": ""},
            {"patient_id": "A2", "age": "60", "sex": "M", "creatinine": "106.08",
             "creatinine_unit": "umol/L", "cystatin_c": "1.0", "acr": "400", "race": "black"},
            {"patient_id": "A3", "age": "70", "sex": "F", "creatinine": "",
             "creatinine_unit": "", "cystatin_c": "", "acr": "", "race": ""},
        ]
        with open(input_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=egfr.CSV_INPUT_FIELDS)
            writer.writeheader()
            writer.writerows(rows)

        results = egfr.process_csv(input_path, output_path)
        assert len(results) == 3

        by_id = {r.patient_id: r for r in results}
        assert by_id["A1"].egfr_2021_creatinine is not None
        assert math.isclose(by_id["A1"].egfr_2021_creatinine, 68.6, abs_tol=0.5)

        assert by_id["A2"].egfr_2021_creatinine is not None
        assert math.isclose(by_id["A2"].egfr_2021_creatinine, 69.2, abs_tol=0.5)
        assert by_id["A2"].albuminuria_stage == "A3"

        assert by_id["A3"].primary_egfr() is None
        assert any("cannot compute" in w for w in by_id["A3"].warnings)

        assert os.path.exists(output_path)
        with open(output_path, encoding="utf-8") as f:
            out_rows = list(csv.DictReader(f))
        assert len(out_rows) == 3


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def test_cli_single():
    rc = egfr.main(["single", "--age", "50", "--sex", "F", "--creatinine", "1.0"])
    assert rc == 0


def test_cli_batch():
    with tempfile.TemporaryDirectory() as tmp:
        inp = os.path.join(tmp, "in.csv")
        out = os.path.join(tmp, "out.csv")
        with open(inp, "w", newline="") as f:
            f.write("patient_id,age,sex,creatinine,creatinine_unit,cystatin_c,acr,race\n")
            f.write("T1,50,F,1.0,mg/dL,,,other\n")
        rc = egfr.main(["batch", "--input", inp, "--output", out])
        assert rc == 0
        assert os.path.exists(out)


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def run_all():
    tests = [obj for name, obj in globals().items() if name.startswith("test_") and callable(obj)]
    passed = 0
    failed = 0
    for t in tests:
        try:
            t()
            passed += 1
            print(f"  PASS: {t.__name__}")
        except Exception as e:
            failed += 1
            print(f"  FAIL: {t.__name__} -- {e}")
    print(f"\n{passed}/{passed + failed} tests passed.")
    return failed


if __name__ == "__main__":
    import sys
    sys.exit(run_all())
