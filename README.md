# Bispecific Antibody Bispecific Binding Agent

> **Domain:** Computational Biology & AI Drug Discovery  
> **Reference Guidelines & Standards:** `wwPDB, IUPAC & CLSI Computational Guidelines`

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

Bispecific Avidity Measurement: avidity coefficient calculation and cooperative binding analysis.

Binding Curve Fitting: global fitting of multi-state binding models with parameter linkage.

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Core Algorithmic & Evaluation Engines

- **`BindingMeasurement`** — dedicated module for binding measurement evaluation and state verification.
- **`AvidityMeasurement`**: Calculate avidity coefficients from bispecific binding data.
- **`BindingCurve`** — dedicated module for binding curve evaluation and state verification.
- **`GlobalBindingFitter`**: Global fitting engine for multi-curve SPR binding data.
- **`TargetEpitope`**: Target antigen epitope specification.
- **`BindingArm`**: Monovalent binding arm kinetics and thermodynamics.

---

## 📐 Mathematical Formulation & Logic

```text
  """Calculate avidity coefficients from bispecific binding data."""
  Calculates optimal antibody concentration [Ab]_opt and maximum ternary complex.
  """Calculates avidity enhancement and effective local search concentration."""
  Calculate effective local concentration C_eff (M) experienced by second arm:
  Calculate apparent avidity KD_app and enhancement factor:
```

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --ab-conc <value> --target-a-conc <value> --target-b-conc <value> --kd-a <value>
```

### Parameter Reference
- `--ab-conc`: Specifies input measurement or parameter value.
- `--target-a-conc`: Specifies input measurement or parameter value.
- `--target-b-conc`: Specifies input measurement or parameter value.
- `--kd-a`: Specifies input measurement or parameter value.
- `--kd-b`: Specifies input measurement or parameter value.
- `--alpha`: Specifies input measurement or parameter value.
- `--json`: Specifies input measurement or parameter value.
- `--linker-len`: Specifies input measurement or parameter value.
- `--c-eff`: Specifies input measurement or parameter value.
- `--construct-id`: Specifies input measurement or parameter value.

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
docker build -t bispecific-antibody-bispecific-binding-agent .
docker run -p 8000:8000 bispecific-antibody-bispecific-binding-agent
```
