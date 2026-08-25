# Bispecific Antibody Bispecific Binding Agent

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python: 3.10+](https://img.shields.io/badge/Python-3.10%2B-brightgreen.svg)](https://python.org)
[![Tests: 30 Passing](https://img.shields.io/badge/Tests-30%20Passing-success.svg)](test_bispecific_antibody.py)
[![Domain: Biophysics & Immuno-oncology](https://img.shields.io/badge/Domain-Biophysics%20%26%20Immuno--Oncology-blueviolet.svg)](#)

A biophysical modeling, kinetics simulation, and rational design platform for **Bispecific Antibodies (BsAbs)**, **Bispecific T-Cell Engagers (BiTEs)**, and **Bivalent Dual-Targeting Biologics**.

---

## Key Biophysical Capabilities

1. **Ternary Complex Mass-Action Equilibrium & Hook Effect**
   - Exact mass-action equilibrium modeling for simultaneous dual-receptor binding:
     $$\text{Ab} + \text{R}_1 \rightleftharpoons \text{Ab}\cdot\text{R}_1 \quad (K_{D1})$$
     $$\text{Ab} + \text{R}_2 \rightleftharpoons \text{Ab}\cdot\text{R}_2 \quad (K_{D2})$$
     $$\text{Ab}\cdot\text{R}_1 + \text{R}_2 \rightleftharpoons \text{Ab}\cdot\text{R}_1\cdot\text{R}_2 \quad (K_{D2}/\alpha)$$
   - Analytical and numeric resolution of bell-shaped ("prozone" / "hook effect") dose-response curves.
   - Calculation of optimal dosing window $[\text{Ab}]_{\text{opt}} \approx \sqrt{K_{D1} K_{D2} (1 + [\text{R}_1]_0/K_{D1})(1 + [\text{R}_2]_0/K_{D2})}$.

2. **Avidity & Effective Local Concentration ($C_{\text{eff}}$)**
   - Polymer physics-based effective search volume calculation:
     $$C_{\text{eff}} = \frac{1}{N_A \cdot V_{\text{search}}} = \frac{3}{4 \pi r_{\text{reach}}^3 N_A}$$
   - Apparent avidity constant $K_{D,\text{app}} = \frac{K_{D1} K_{D2}}{C_{\text{eff}} + K_{D2}}$ and avidity enhancement ratio $\beta = \frac{K_{D1}}{K_{D,\text{app}}}$.

3. **BiTE / T-Cell Engager Cytotoxicity & Synapse Modeling**
   - Effector-to-Target (E:T) ratio and immunological synapse formation kinetics.
   - Specific target cell lysis percentage predictions via Hill cooperativity sigmoidal dose-response models.
   - Off-target cytokine release syndrome (CRS) risk index quantification.

4. **Multi-Cycle Surface Plasmon Resonance (SPR) Global Regression**
   - 1:1 Langmuir binding sensorgram simulations (association & dissociation phases).
   - Log-scale parameterization non-linear multi-concentration global fitting of $k_a$ ($\text{M}^{-1}\text{s}^{-1}$), $k_d$ ($\text{s}^{-1}$), $K_D$, and $R_{\max}$ with $R^2 > 0.98$ and $\chi^2$ goodness-of-fit validation.

5. **Worm-Like Chain (WLC) Linker Mechanics Optimization**
   - Flexible peptide linker $((Gly_4Ser)_n)$ modeling via Kratky-Porod polymer mechanics:
     $$\langle R^2 \rangle^{1/2} = \sqrt{2 l_p L_c \left(1 - \frac{l_p}{L_c}(1 - e^{-L_c / l_p})\right)}$$
   - Geometric strain energy and structural reach feasibility analysis across candidate epitope distances.

---

## Installation

```bash
git clone https://github.com/example/bispecific-antibody-bispecific-binding-agent.git
cd bispecific-antibody-bispecific-binding-agent
```

*Requires Python 3.10+ with zero external third-party dependencies (pure standard library).*

---

## Command-Line Interface (CLI)

```bash
# 1. Solve ternary complex equilibrium for a bispecific construct
python cli.py ternary --ab-conc 1e-9 --target-a-conc 1e-8 --target-b-conc 1e-8 --kd-a 1e-9 --kd-b 1e-8

# 2. Compute dose-response titration curve and detect hook effect
python cli.py titration --target-a-conc 1e-8 --target-b-conc 1e-8 --kd-a 1e-9 --kd-b 1e-8

# 3. Calculate bivalent avidity enhancement & effective concentration C_eff
python cli.py avidity --kd-a 1e-8 --kd-b 1e-8 --linker-len 15

# 4. Predict T-cell mediated target lysis and CRS risk
python cli.py tcell-lysis --tumor-density 50000 --cd3-density 40000 --et-ratio 5.0 --ab-conc 1e-10

# 5. Run global multi-concentration SPR curve fitting
python cli.py fit-spr

# 6. Optimize peptide linker length for target epitope distance
python cli.py optimize-linker --distance 4.5 --kd-a 1e-9 --kd-b 1e-8

# 7. Output in JSON format (available on all subcommands)
python cli.py ternary --ab-conc 1e-9 --json
```

---

## Python API Usage

```python
from bispecific_antibody import (
    TernaryComplexModel,
    AvidityEngine,
    TCellEngagerPredictor,
    BispecificConstruct,
    BindingArm,
    LinkerProperties,
    LinkerOptimizer,
)

# Solve ternary complex formation at equilibrium
eq = TernaryComplexModel.solve_equilibrium(
    antibody_conc_m=1e-9,
    target_a_conc_m=1e-8,
    target_b_conc_m=1e-8,
    kd_a=1e-9,
    kd_b=1e-8,
    alpha=1.0,
)
print(f"Ternary Complex: {eq['ternary_complex']:.4e} M")

# Calculate bivalent avidity enhancement
avidity = AvidityEngine.calculate_apparent_avidity(kd_a=1e-8, kd_b=1e-8, linker_length_aa=15)
print(f"Apparent KD: {avidity['apparent_avidity_kd_m']:.4e} M (Enhancement: {avidity['avidity_enhancement_factor']:.2f}x)")

# Simulate BiTE-mediated cytotoxicity
construct = BispecificConstruct(
    construct_id="BiTE-CD19xCD3",
    name="Blinatumomab-like",
    arm_a=BindingArm("CD19", ka=1e5, kd=1e-4),
    arm_b=BindingArm("CD3", ka=1e5, kd=1e-3),
    linker=LinkerProperties("(G4S)3", length_aa=15),
)
lysis = TCellEngagerPredictor.predict_lysis(construct, target_cell_density=50000.0, antibody_conc_m=1e-10)
print(f"Tumor Lysis: {lysis['specific_lysis_percentage']:.1f}% | CRS Risk: {lysis['crs_risk_category']}")
```

---

## Test Suite

Execute the unit test suite:

```bash
python -m unittest test_bispecific_antibody.py
```

All 30 unit tests cover:
- Zero, asymmetric, and boundary concentration edge cases
- Positive/negative cooperativity ($\alpha$) scaling
- High-dose hook effect attenuation
- WLC polymer mechanics and RMS span
- Effector-to-target ratio and cytotoxicity sigmoids
- Global SPR sensorgram fitting regression
- CLI subcommand JSON and human-readable output validation

---

## License

MIT License. See `LICENSE` for details.
