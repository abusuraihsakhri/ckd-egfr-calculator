# CKD-EPI eGFR Calculator

Estimates glomerular filtration rate (eGFR) using the CKD-EPI equations and stages chronic kidney disease per KDIGO 2012 guidelines.

## Equations Implemented

### CKD-EPI 2021 (Race-Free)
- **Creatinine**: `eGFR = 142 × min(Scr/κ, 1)^α × max(Scr/κ, 1)^(-1.200) × 0.9938^Age × sex_factor`
  - κ = 0.7 (F) / 0.9 (M), α = -0.241 (F) / -0.302 (M), sex_factor = 1.012 (F) / 1.0 (M)
- **Cystatin C**: `eGFR = 133 × min(Scys/0.8, 1)^(-0.499) × max(Scys/0.8, 1)^(-1.328) × 0.9961^Age × sex_factor`
- **Combined (Cr-CysC)**: Uses both markers for highest accuracy

### CKD-EPI 2009 (With Race Coefficient)
- Includes the 1.159 race multiplier for Black patients (Levey et al., Ann Intern Med 2009)
- Provided for reference/comparison purposes

### KDIGO 2012 Staging
- **GFR Categories**: G1 (≥90), G2 (60-89), G3a (45-59), G3b (30-44), G4 (15-29), G5 (<15)
- **Albuminuria Categories**: A1 (<30), A2 (30-300), A3 (>300 mg/g)
- **Combined Risk Grid**: Low → Very High risk based on G + A stage

## Usage

```bash
# Single patient
python egfr_calculator.py single --age 55 --sex M --creatinine 1.2

# With all markers
python egfr_calculator.py single --age 65 --sex F --creatinine 1.3 --cystatin-c 1.1 --acr 150

# With race for 2009 equation
python egfr_calculator.py single --age 50 --sex M --creatinine 1.0 --race black

# Batch CSV processing
python egfr_calculator.py batch --input patients.csv --output results.csv
```

## CSV Input Format

Required columns: `patient_id`, `age`, `sex`
Optional columns: `creatinine`, `creatinine_unit` (mg/dL or umol/L), `cystatin_c`, `acr`, `race`

## Requirements

Python 3.9+ (stdlib only, no external dependencies)

## Disclaimer

This calculator is for educational and clinical decision support purposes only. It does not replace professional medical judgment. The 2021 CKD-EPI equations are race-free per current NKF-ASN guidance. The 2009 equation with race coefficient is included for historical reference.
