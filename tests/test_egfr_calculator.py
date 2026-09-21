import csv
import math
import os
import tempfile

import pytest
import egfr_calculator as egfr


def test_2021_creatinine_reference_case():
    result = egfr.ckd_epi_2021_creatinine(60, "M", 1.2)
    assert math.isclose(result, 69.2, abs_tol=0.5)


def test_2012_cystatin_formula_uses_2012_age_factor():
    expected = 133 * 0.996**50 * 0.932
    assert math.isclose(egfr.ckd_epi_2012_cystatin(50, "F", 0.8), expected, abs_tol=1e-12)
    assert egfr.ckd_epi_2021_cystatin(50, "F", 0.8) == egfr.ckd_epi_2012_cystatin(50, "F", 0.8)


def test_2021_combined_reference_formula():
    result = egfr.ckd_epi_2021_creatinine_cystatin(65, "F", 1.3, 1.1)
    expected = 135 * (1.3 / 0.7) ** -0.544 * (1.1 / 0.8) ** -0.778 * 0.9961**65 * 0.963
    assert math.isclose(result, expected, abs_tol=1e-12)


@pytest.mark.parametrize(
    "value,expected",
    [(90, "G1"), (89.9999995, "G2"), (60, "G2"), (59.9999995, "G3a"),
     (45, "G3a"), (44.9999995, "G3b"), (30, "G3b"), (29.9999995, "G4"),
     (15, "G4"), (14.9999995, "G5")],
)
def test_gfr_boundaries_have_no_floating_point_gaps(value, expected):
    assert egfr.classify_gfr_stage(value)[0] == expected


@pytest.mark.parametrize("value,expected", [(29.9999999, "A1"), (30, "A2"), (299.9999999, "A2"), (300, "A3")])
def test_albuminuria_boundaries(value, expected):
    assert egfr.classify_albuminuria(value)[0] == expected


@pytest.mark.parametrize("func,args", [
    (egfr.classify_gfr_stage, (-1,)),
    (egfr.classify_gfr_stage, (float("nan"),)),
    (egfr.classify_albuminuria, (-1,)),
    (egfr.ckd_epi_2021_creatinine, (17, "M", 1.0)),
    (egfr.ckd_epi_2021_creatinine, (50, "M", 0)),
    (egfr.ckd_epi_2012_cystatin, (50, "F", -1)),
])
def test_invalid_inputs_reject(func, args):
    with pytest.raises(ValueError):
        func(*args)


def test_unit_conversion_accepts_micro_symbol():
    assert math.isclose(egfr.creatinine_to_mgdl(88.4, "µmol/L"), 1.0)
    assert math.isclose(egfr.creatinine_to_mgdl(88.4, "umol/L"), 1.0)


def test_patient_combined_preferred_and_categories():
    result = egfr.calculate_patient("P", 65, "F", 1.3, 1.1, 150)
    assert result.primary_egfr() == result.egfr_2021_combined
    assert result.egfr_2012_cystatin is not None
    assert result.egfr_2021_cystatin == result.egfr_2012_cystatin
    assert result.albuminuria_stage == "A2"
    assert result.combined_risk is not None


def test_batch_bad_optional_value_does_not_abort_other_rows():
    with tempfile.TemporaryDirectory() as tmp:
        inp = os.path.join(tmp, "input.csv")
        out = os.path.join(tmp, "output.csv")
        with open(inp, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(egfr.CSV_INPUT_FIELDS)
            writer.writerow(["P1", "50", "F", "1.0", "mg/dL", "bad", "20", "other"])
            writer.writerow(["P2", "60", "M", "1.2", "mg/dL", "1.0", "400", "other"])
        results = egfr.process_csv(inp, out)
        assert len(results) == 2
        assert results[0].egfr_2021_creatinine is not None
        assert any("cystatin C" in warning for warning in results[0].warnings)
        assert results[1].egfr_2021_combined is not None
        with open(out, encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        assert "egfr_2012_cystatin" in rows[0]
        assert "gfr_category" in rows[0]


def test_cli_returns_error_for_pediatric_input():
    assert egfr.main(["single", "--age", "12", "--sex", "M", "--creatinine", "1.0"]) == 2
