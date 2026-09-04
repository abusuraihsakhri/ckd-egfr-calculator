# CKD-EPI eGFR Calculator & KDIGO Staging Engine

> **Domain:** Nephrology, Renal Epidemiology & Clinical Laboratory Diagnostics  
> **Clinical Guidelines & Standards:** 2021 CKD-EPI Race-Free Equations (Inker et al., NEJM 2021), 2009 CKD-EPI Equation (Levey et al., Ann Intern Med 2009), KDIGO 2012 / 2024 Clinical Practice Guideline for the Evaluation and Management of Chronic Kidney Disease

---

## 📖 Clinical Overview

The **CKD-EPI eGFR Calculator** computes estimated glomerular filtration rate (eGFR) and risk-stratifies chronic kidney disease (CKD) using the validated KDIGO 2012 / 2024 heat map. It implements both the updated 2021 race-free CKD-EPI equations (Creatinine, Cystatin C, and Combined Cr-CysC) and the legacy 2009 CKD-EPI equation for longitudinal trend comparison.

### Mathematical Formulations

#### 1. 2021 Race-Free CKD-EPI Creatinine Equation
$$\text{eGFR}_{\text{cr}} = 142 \times \min\left(\frac{S_{\text{cr}}}{\kappa}, 1\right)^\alpha \times \max\left(\frac{S_{\text{cr}}}{\kappa}, 1\right)^{-1.200} \times 0.9938^{\text{Age}} \times [1.012 \text{ if Female}]$$
Where:
- Female: $\kappa = 0.7$, $\alpha = -0.241$
- Male: $\kappa = 0.9$, $\alpha = -0.302$

#### 2. 2021 Race-Free CKD-EPI Creatinine-Cystatin C Equation
$$\begin{aligned}
\text{eGFR}_{\text{cr-cys}} = & 135 \times \min\left(\frac{S_{\text{cr}}}{\kappa}, 1\right)^\alpha \times \max\left(\frac{S_{\text{cr}}}{\kappa}, 1\right)^{-0.544} \times \min\left(\frac{S_{\text{cys}}}{0.8}, 1\right)^{-0.323} \\
& \times \max\left(\frac{S_{\text{cys}}}{0.8}, 1\right)^{-0.778} \times 0.9961^{\text{Age}} \times [0.963 \text{ if Female}]
\end{aligned}$$
Where:
- Female: $\kappa = 0.7$, $\alpha = -0.219$
- Male: $\kappa = 0.9$, $\alpha = -0.144$

### KDIGO CKD Staging Framework
- **GFR Categories:** G1 ($\ge 90$), G2 (60–89), G3a (45–59), G3b (30–44), G4 (15–29), G5 ($< 15\,\text{mL/min/1.73m}^2$, Kidney Failure).
- **Albuminuria Categories (ACR in mg/g):** A1 ($< 30$, Normal/Mild), A2 (30–300, Moderate), A3 ($> 300$, Severe).

---

## 💻 CLI Quickstart & Usage

### 1. Calculate Single Patient eGFR & Staging
```bash
python cli.py single --age 60 --sex M --creatinine 1.2 --cystatin-c 1.0 --acr 400
```

### 2. Batch Process Patient CSV Dataset
```bash
python cli.py batch -i sample.csv -o out_results.csv
```

---

## 🧪 Verification & Testing

Execute comprehensive unit tests via pytest:
```bash
python -m pytest -p no:zarr
```
