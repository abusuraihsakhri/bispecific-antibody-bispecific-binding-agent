# Bispecific Antibody Bispecific Binding Agent & Kinetics Engine

> **Domain:** Computational Biophysics, Quantitative Pharmacology & AI Drug Discovery  
> **Reference Guidelines & Standards:** `wwPDB, IUPAC & CLSI Computational Guidelines, MIQE/SPR Standards`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![Pytest](https://img.shields.io/badge/Pytest-Passing-brightgreen.svg)
![Biophysics](https://img.shields.io/badge/Model-Ternary_Equilibrium-orange.svg)

</div>

---

## 📖 Overview & Scientific Scope

The **Bispecific Antibody Kinetics Engine** provides rigorous mathematical and biophysical simulations of bispecific antibody constructs (including T-cell engagers/BiTEs, dual tumor-targeting IgG antibodies, and cross-linking scaffolds).

The platform quantifies:
1. **Ternary Complex Equilibrium $[A \cdot \text{Ab} \cdot B]$**: Non-linear mass action equilibrium between bispecific antibodies ($\text{Ab}$) and two distinct targets or cellular receptors ($A$ and $B$).
2. **Bell-Shaped Hook Effect (Prozone Phenomenon)**: Biphasic dose-response curve wherein supra-optimal antibody concentrations favor binary complexes over functional bridging ternary complexes.
3. **Avidity Enhancement & Effective Local Concentration ($C_{\text{eff}}$)**: Thermodynamics of bivalent dual-epitope binding governed by flexible peptide linkers via polymer physics (Worm-Like Chain model).
4. **T-Cell Immunological Synapse Formation & Cytotoxicity Prediction**: Receptor cross-linking densities, target cell lysis percentage, and Cytokine Release Syndrome (CRS) risk grading.
5. **Global Multi-Cycle SPR Kinetics Fitting**: Global Langmuir 1:1 regression across multi-concentration surface plasmon resonance sensorgrams.

---

## 📐 Mathematical Formulation & Biophysical Modeling

### 1. Ternary Complex Equilibrium $[A \cdot \text{Ab} \cdot B]$

For a bispecific antibody $\text{Ab}$ binding two antigen targets $A$ and $B$, the system obeys the following coupled equilibria:

$$\text{Ab} + A \xrightleftharpoons{K_{D,A}} [\text{Ab} \cdot A], \quad K_{D,A} = \frac{[\text{Ab}]_{\text{free}} [A]_{\text{free}}}{[\text{Ab} \cdot A]}$$

$$\text{Ab} + B \xrightleftharpoons{K_{D,B}} [\text{Ab} \cdot B], \quad K_{D,B} = \frac{[\text{Ab}]_{\text{free}} [B]_{\text{free}}}{[\text{Ab} \cdot B]}$$

$$[\text{Ab} \cdot A] + B \xrightleftharpoons{K_{D,B}/\alpha} [A \cdot \text{Ab} \cdot B], \quad [\text{Ab} \cdot B] + A \xrightleftharpoons{K_{D,A}/\alpha} [A \cdot \text{Ab} \cdot B]$$

where $\alpha$ is the thermodynamic **cooperativity factor**:
- $\alpha = 1$: Independent, non-cooperative binding.
- $\alpha > 1$: Positive allosteric/avidity cooperativity (enhanced ternary formation).
- $\alpha < 1$: Negative cooperativity (steric clash or antagonistic conformational change).

The trimeric complex concentration $[A \cdot \text{Ab} \cdot B]$ is given by:

$$[A \cdot \text{Ab} \cdot B] = \frac{[\text{Ab}]_{\text{free}} [A]_{\text{free}} [B]_{\text{free}}}{K_{D,A} (K_{D,B} / \alpha)}$$

### 2. The Hook Effect (Biphasic Titration Curve)

At low antibody concentrations, ternary complex formation increases monotonically with $[\text{Ab}]_0$. However, once antibody concentration exceeds target receptor saturation, competitive binary binding takes over, yielding two binary complexes ($[\text{Ab} \cdot A]$ and $[\text{Ab} \cdot B]$) rather than one ternary bridge:

$$[\text{Ab}]_{\text{opt}} \approx \sqrt{K_{D,A} \cdot K_{D,B} \left(1 + \frac{[A]_0}{K_{D,A}}\right) \left(1 + \frac{[B]_0}{K_{D,B}}\right)}$$

Supra-optimal dosing leads to significant loss of therapeutic efficacy (the classical immunological hook / prozone effect).

### 3. Dual Avidity Enhancement & Effective Local Concentration ($C_{\text{eff}}$)

When an antibody binds cell-surface epitope $A$ with arm 1, arm 2 experiences an effective local concentration $C_{\text{eff}}$ confined within the search volume of the flexible linker:

$$V_{\text{search}} = \frac{4}{3} \pi \langle R^2 \rangle^{3/2}$$

$$C_{\text{eff}} = \frac{1}{N_A \cdot V_{\text{search}}}$$

where $\langle R^2 \rangle^{1/2}$ is the Root-Mean-Square (RMS) end-to-end distance calculated from the **Worm-Like Chain (WLC)** model:

$$\langle R^2 \rangle = 2 l_p L_c \left[ 1 - \frac{l_p}{L_c} \left( 1 - e^{-L_c / l_p} \right) \right]$$

with contour length $L_c = N_{\text{aa}} \times 0.38\text{ nm}$ and peptide persistence length $l_p \approx 0.4\text{ nm}$.

The apparent bivalent avidity $K_{D,\text{app}}$ and enhancement factor $\beta$ are:

$$K_{D,\text{app}} = \frac{K_{D,A} \cdot K_{D,B}}{C_{\text{eff}} + K_{D,B}}$$

$$\beta = \frac{K_{D,A}}{K_{D,\text{app}}} = 1 + \frac{C_{\text{eff}}}{K_{D,B}}$$

---

## 💻 CLI Quickstart & Usage

The CLI tool `cli.py` provides high-performance simulation subcommands.

### 1. Batch Screening from CSV (`sample.csv`)
Process multiple bispecific antibody designs, T-cell engagers, and dual-targeting constructs:

```bash
# Run batch analysis and export results to output CSV
python cli.py batch -i sample.csv -o batch_results.csv

# Output results as JSON directly to stdout
python cli.py batch -i sample.csv --json
```

### 2. Ternary Complex Equilibrium
Calculate species concentrations at specific antibody and target receptor concentrations:

```bash
python cli.py ternary \
  --ab-conc 1e-10 \
  --target-a-conc 5e-9 \
  --target-b-conc 5e-9 \
  --kd-a 1.5e-9 \
  --kd-b 1e-7 \
  --alpha 1.0
```

### 3. Hook Effect Titration Profile
Evaluate the dose-response curve and detect bell-shaped attenuation:

```bash
python cli.py titration \
  --target-a-conc 1e-8 \
  --target-b-conc 1e-8 \
  --kd-a 1e-9 \
  --kd-b 1e-8 \
  --alpha 1.0
```

### 4. Bivalent Avidity & Linker Evaluation
Calculate effective local concentration $C_{\text{eff}}$ and apparent avidity $K_{D,\text{app}}$:

```bash
python cli.py avidity \
  --kd-a 1e-8 \
  --kd-b 1e-8 \
  --linker-len 15
```

### 5. T-Cell Engager (BiTE) Cytotoxicity Simulation
Simulate immunological synapse density per tumor cell, specific lysis percentage, and CRS risk:

```bash
python cli.py tcell-lysis \
  --construct-id BiTE-CD19 \
  --tumor-antigen CD19 \
  --tumor-kd 1.5e-9 \
  --cd3-kd 1.0e-7 \
  --tumor-density 50000 \
  --cd3-density 40000 \
  --et-ratio 5.0 \
  --ab-conc 1e-10
```

### 6. Polymer Physics Linker Optimization
Identify the optimal linker length bridging two target epitopes at known spatial distance:

```bash
python cli.py optimize-linker \
  --distance 4.5 \
  --kd-a 1e-9 \
  --kd-b 1e-8
```

---

## 🐍 Python API Quickstart

```python
from bispecific_antibody import (
    TernaryComplexModel,
    AvidityEngine,
    BispecificConstruct,
    BindingArm,
    LinkerProperties,
    TCellEngagerPredictor,
)

# 1. Calculate ternary complex equilibrium
ternary = TernaryComplexModel.solve_equilibrium(
    antibody_conc_m=1e-10,
    target_a_conc_m=5e-9,
    target_b_conc_m=5e-9,
    kd_a=1.5e-9,
    kd_b=1e-7,
    alpha=1.0,
)
print("Ternary complex concentration:", ternary["ternary_complex"])
print("Target A occupancy:", ternary["ternary_fraction_of_target_a"])

# 2. Avidity enhancement calculation
avidity = AvidityEngine.calculate_apparent_avidity(
    kd_a=1.5e-9,
    kd_b=1e-7,
    linker_length_aa=15,
)
print("Effective Concentration C_eff:", avidity["effective_concentration_m"])
print("Apparent Avidity KD_app:", avidity["apparent_avidity_kd_m"])
print("Enhancement Factor:", avidity["avidity_enhancement_factor"])

# 3. T-cell cytotoxicity simulation
construct = BispecificConstruct(
    construct_id="BSAB-001",
    name="Blinatumomab-CD19xCD3",
    arm_a=BindingArm("CD19", ka=1e5, kd=1.5e-4),
    arm_b=BindingArm("CD3", ka=1e5, kd=1.0e-2),
    linker=LinkerProperties(sequence_or_type="(G4S)3", length_aa=15),
)
lysis = TCellEngagerPredictor.predict_lysis(
    construct=construct,
    target_cell_density=50000,
    cd3_density=40000,
    effector_to_target_ratio=5.0,
    antibody_conc_m=1e-10,
)
print("Synapses / cell:", lysis["synapses_per_tumor_cell"])
print("Specific Lysis %:", lysis["specific_lysis_percentage"])
print("CRS Risk:", lysis["crs_risk_category"])
```

---

## 🧪 Testing & Verification

Run the full automated test suite:

```bash
python -m pytest -p no:zarr -v
```

Verify CLI batch execution:

```bash
python cli.py batch -i sample.csv -o out_smoke.csv
```
