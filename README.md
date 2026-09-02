# CKD EGFR Calculator

> **Domain:** Nephrology & Renal Replacement Protocols  
> **Reference Guidelines & Standards:** `KDIGO & KDOQI Clinical Guidelines`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

CKD-EPI eGFR Calculator
========================

Computes estimated glomerular filtration rate (eGFR) using:
  - 2021 CKD-EPI creatinine equation (race-free, Inker et al., NEJM 2021)
  - 2021 CKD-EPI cystatin C equation (race-free)
  - 2021 CKD-EPI creatinine-cystatin C combined equation (race-free)
  - 2009 CKD-EPI creatinine equation (with race coefficient, Levey et al.)

Stages chronic kidney disease per KDIGO 2012 GFR/albuminuria grid.

Stdlib only. Usage: python egfr_calculator.py --help

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Core Algorithmic & Evaluation Engines

- **`EgfrResult`** — dedicated module for egfr result evaluation and state verification.

---

## 📐 Mathematical Formulation & Logic

```text
  """Return (stage_label, description) for a KDIGO GFR category."""
  """Return (stage_label, description) for a KDIGO albuminuria category."""
  result.combined_risk = kdigo_combined_risk(stage, a_stage)
  result = calculate_patient(
  single = subparsers.add_parser("single", help="Calculate eGFR for one patient")
```

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --input data.csv
```

### Parameter Reference
- `--interactive`: Launch guided terminal interactive wizard.
- `--input <path>`: Evaluate input from JSON or CSV specification.
- `--json`: Output deterministic structured results in JSON format.

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `suite_name` | Parameter / observation metric | Required |
| `system_slug` | Parameter / observation metric | Required |
| `standard_reference` | Parameter / observation metric | Required |
| `test_cases` | Parameter / observation metric | Required |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t ckd-egfr-calculator .
docker run -p 8000:8000 ckd-egfr-calculator
```
