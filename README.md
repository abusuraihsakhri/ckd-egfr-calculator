# CKD-EPI eGFR Calculator

A compact adult eGFR calculator implementing the 2021 CKD-EPI creatinine and creatinine–cystatin C equations, the 2012 CKD-EPI cystatin C equation, and KDIGO GFR/albuminuria categories. The repository also provides a dependency-free Python CLI and CSV batch workflow.

## Features

- 2021 race-free CKD-EPI creatinine equation (`eGFRcr`)
- 2012 CKD-EPI cystatin C equation (`eGFRcys`)
- 2021 race-free CKD-EPI creatinine–cystatin C equation (`eGFRcr-cys`)
- Legacy 2009 CKD-EPI creatinine calculation in the Python CLI for longitudinal comparison only
- KDIGO G1–G5 and A1–A3 categorization with the G/A risk grid
- Creatinine input in mg/dL or µmol/L
- Browser interface with light/dark themes and no server-side processing
- CSV batch processing from the command line

## Use

Open the live application, enter age and sex, then provide creatinine, cystatin C, or both. ACR is optional and enables albuminuria categorization and the combined KDIGO risk category. Press **Analyze** to calculate.

For the CLI:

```bash
python cli.py single --age 60 --sex M --creatinine 1.2 --cystatin-c 1.0 --acr 400
python cli.py batch -i sample.csv -o results.csv
```

The input CSV uses these columns:

```text
patient_id,age,sex,creatinine,creatinine_unit,cystatin_c,acr,race
```

`race` is used only for the retained 2009 legacy comparison. Current primary estimates use race-free equations.

## Clinical interpretation

This tool reports estimated GFR and KDIGO GFR/albuminuria categories. A single measurement does **not** diagnose chronic kidney disease. CKD requires evidence of chronicity and/or other markers of kidney damage, and results should be interpreted with the clinical context and appropriately standardized assays.

The formulas implemented here correspond to the National Kidney Foundation/NIDDK descriptions of the 2021 CKD-EPI creatinine and creatinine–cystatin C equations and the 2012 CKD-EPI cystatin C equation. KDIGO 2024 retains the G1–G5 and A1–A3 classification framework.

## Privacy

The web calculator runs entirely in the browser. Patient values are not sent to a server or persisted by the application. Only the selected light/dark theme may be stored locally. The Python CLI processes files locally.

## Development and testing

Runtime code uses only the Python standard library and browser-native JavaScript. `pytest` is needed only for Python tests.

```bash
python -m pip install pytest
python -m pytest -q
node --test tests/web_calculator.test.mjs
python -m compileall -q egfr_calculator.py cli.py
```

GitHub Actions runs Python tests on Python 3.10–3.13, browser-core tests on Node 22, a CLI smoke test, and GitHub Pages deployment.

## Browser compatibility

The web app uses standard ES modules and modern CSS. Current versions of Chrome, Edge, Firefox, and Safari are supported. JavaScript must be enabled.

## License

MIT License. See [LICENSE](LICENSE).
