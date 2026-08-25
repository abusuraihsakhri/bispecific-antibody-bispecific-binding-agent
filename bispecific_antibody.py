"""
Bispecific Antibody Binding Agent & Kinetics Engine
====================================================
Comprehensive biophysical modeling for bispecific antibodies (BsAbs, BiTEs, Dual-IgGs):
- Trimeric/ternary complex equilibrium modeling with bell-shaped hook effect
- Avidity and effective concentration (C_eff) thermodynamics
- T-cell immunological synapse formation and cytotoxicity dose-response
- Global SPR kinetics simulation and non-linear multi-concentration curve fitting
- Polymer physics-based linker mechanics (Worm-Like Chain model)
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple


@dataclass
class TargetEpitope:
    """Target antigen epitope specification."""
    name: str
    cell_type: str  # e.g., "Tumor", "T-Cell", "NK-Cell"
    density_per_cell: float  # molecules per cell
    target_concentration_m: float  # molar concentration in compartment (M)


@dataclass
class BindingArm:
    """Monovalent binding arm kinetics and thermodynamics."""
    target_name: str
    ka: float  # on-rate (M^-1 s^-1)
    kd: float  # off-rate (s^-1)
    valency: int = 1

    @property
    def kd_molar(self) -> float:
        """Equilibrium dissociation constant KD (M)."""
        return self.kd / self.ka if self.ka > 0 else float("inf")


@dataclass
class LinkerProperties:
    """Linker structural parameters."""
    sequence_or_type: str  # e.g. "(G4S)3"
    length_aa: int
    persistence_length_nm: float = 0.4  # ~0.4 nm for flexible peptide (Gly-Ser)
    aa_step_nm: float = 0.38  # 0.38 nm per peptide bond in extended state

    @property
    def contour_length_nm(self) -> float:
        """Total contour length L_c = N_aa * 0.38 nm."""
        return self.length_aa * self.aa_step_nm

    @property
    def root_mean_square_end_to_end_nm(self) -> float:
        """Worm-Like Chain (Kratky-Porod) RMS end-to-end distance <R^2>^0.5."""
        lc = self.contour_length_nm
        lp = self.persistence_length_nm
        if lc <= 0 or lp <= 0:
            return 0.0
        ratio = lp / lc
        r2 = 2.0 * lp * lc * (1.0 - ratio * (1.0 - math.exp(-lc / lp)))
        return math.sqrt(max(0.0, r2))


@dataclass
class BispecificConstruct:
    """Bispecific antibody construct."""
    construct_id: str
    name: str
    arm_a: BindingArm
    arm_b: BindingArm
    linker: LinkerProperties
    cooperativity_alpha: float = 1.0  # cooperativity factor alpha (>1 positive, <1 negative)
    molecular_weight_kda: float = 55.0  # e.g., 55 kDa for BiTE, 150 kDa for full BsAb


@dataclass
class SPRCurve:
    """SPR sensorgram dataset for a given analyte concentration."""
    concentration_m: float
    time_s: List[float]
    response_ru: List[float]
    association_duration_s: float


# ==============================================================================
# 1. TERNARY COMPLEX EQUILIBRIUM & HOOK EFFECT
# ==============================================================================

class TernaryComplexModel:
    """
    Mathematical model for bispecific antibody ternary complex formation.
    Reactions:
      Ab + TargetA <-> Ab:TargetA (KD_A)
      Ab + TargetB <-> Ab:TargetB (KD_B)
      Ab:TargetA + TargetB <-> TargetA:Ab:TargetB (KD_B / alpha)
      Ab:TargetB + TargetA <-> TargetA:Ab:TargetB (KD_A / alpha)
    """

    @staticmethod
    def solve_equilibrium(
        antibody_conc_m: float,
        target_a_conc_m: float,
        target_b_conc_m: float,
        kd_a: float,
        kd_b: float,
        alpha: float = 1.0,
        max_iter: int = 100,
        tolerance: float = 1e-12,
    ) -> Dict[str, float]:
        """
        Solve equilibrium mass-action equations for free and bound species.
        Returns concentrations of:
          - free_ab, free_a, free_b
          - binary_ab_a, binary_ab_b
          - ternary_complex
        """
        if antibody_conc_m < 0 or target_a_conc_m < 0 or target_b_conc_m < 0:
            raise ValueError("Concentrations must be non-negative")
        if kd_a <= 0 or kd_b <= 0 or alpha <= 0:
            raise ValueError("KD and cooperativity alpha must be positive")

        if antibody_conc_m == 0 or target_a_conc_m == 0 or target_b_conc_m == 0:
            return {
                "free_ab": antibody_conc_m,
                "free_a": target_a_conc_m,
                "free_b": target_b_conc_m,
                "binary_ab_a": 0.0,
                "binary_ab_b": 0.0,
                "ternary_complex": 0.0,
                "ternary_fraction_of_target_a": 0.0,
                "ternary_fraction_of_target_b": 0.0,
            }

        kd_a_eff = kd_a / alpha
        kd_b_eff = kd_b / alpha

        # Damped fixed-point iteration for free receptor and antibody concentrations
        free_a = target_a_conc_m
        free_b = target_b_conc_m
        free_ab = antibody_conc_m

        for _ in range(max_iter):
            ab_denom = 1.0 + (free_a / kd_a) + (free_b / kd_b) + (free_a * free_b / (kd_a * kd_b_eff))
            new_free_ab = antibody_conc_m / max(ab_denom, 1e-30)

            a_denom = 1.0 + (new_free_ab / kd_a) + (new_free_ab * free_b / (kd_a * kd_b_eff))
            new_free_a = target_a_conc_m / max(a_denom, 1e-30)

            b_denom = 1.0 + (new_free_ab / kd_b) + (new_free_ab * new_free_a / (kd_a * kd_b_eff))
            new_free_b = target_b_conc_m / max(b_denom, 1e-30)

            delta = abs(new_free_ab - free_ab) + abs(new_free_a - free_a) + abs(new_free_b - free_b)

            free_ab = 0.5 * free_ab + 0.5 * new_free_ab
            free_a = 0.5 * free_a + 0.5 * new_free_a
            free_b = 0.5 * free_b + 0.5 * new_free_b

            if delta < tolerance:
                break

        binary_ab_a = (free_ab * free_a) / kd_a
        binary_ab_b = (free_ab * free_b) / kd_b
        ternary = (free_ab * free_a * free_b) / (kd_a * kd_b_eff)

        frac_a = ternary / target_a_conc_m if target_a_conc_m > 0 else 0.0
        frac_b = ternary / target_b_conc_m if target_b_conc_m > 0 else 0.0

        return {
            "free_ab": free_ab,
            "free_a": free_a,
            "free_b": free_b,
            "binary_ab_a": binary_ab_a,
            "binary_ab_b": binary_ab_b,
            "ternary_complex": ternary,
            "ternary_fraction_of_target_a": frac_a,
            "ternary_fraction_of_target_b": frac_b,
        }

    @classmethod
    def compute_dose_response_curve(
        cls,
        target_a_conc_m: float,
        target_b_conc_m: float,
        kd_a: float,
        kd_b: float,
        alpha: float = 1.0,
        ab_conc_range: Optional[List[float]] = None,
    ) -> Dict[str, Any]:
        """
        Generate a multi-point titration curve showing bell-shaped hook effect.
        Calculates optimal antibody concentration [Ab]_opt and maximum ternary complex.
        """
        if ab_conc_range is None:
            geom_kd = math.sqrt(kd_a * kd_b)
            ab_conc_range = [geom_kd * (10 ** (x / 4.0)) for x in range(-24, 25)]

        points = []
        max_ternary = -1.0
        optimal_ab_conc = 0.0

        for ab_conc in ab_conc_range:
            res = cls.solve_equilibrium(ab_conc, target_a_conc_m, target_b_conc_m, kd_a, kd_b, alpha)
            ternary = res["ternary_complex"]
            if ternary > max_ternary:
                max_ternary = ternary
                optimal_ab_conc = ab_conc
            points.append({
                "antibody_conc_m": ab_conc,
                "ternary_complex_m": ternary,
                "binary_ab_a_m": res["binary_ab_a"],
                "binary_ab_b_m": res["binary_ab_b"],
                "free_ab_m": res["free_ab"],
            })

        theoretical_opt = math.sqrt(kd_a * kd_b * (1.0 + target_a_conc_m / kd_a) * (1.0 + target_b_conc_m / kd_b))

        return {
            "dose_response_points": points,
            "max_ternary_complex_m": max_ternary,
            "optimal_antibody_conc_m": optimal_ab_conc,
            "theoretical_optimal_ab_m": theoretical_opt,
            "hook_effect_present": len(points) > 2 and points[-1]["ternary_complex_m"] < max_ternary * 0.5,
        }


# ==============================================================================
# 2. AVIDITY & EFFECTIVE LOCAL CONCENTRATION
# ==============================================================================

class AvidityEngine:
    """Calculates avidity enhancement and effective local search concentration."""

    AVOGADRO = 6.02214076e23

    @classmethod
    def calculate_effective_concentration(cls, span_nm: float) -> float:
        """
        Calculate effective local concentration C_eff (M) experienced by second arm:
        C_eff = 1 / (N_A * V_sphere) where V_sphere = 4/3 * pi * r^3
        """
        if span_nm <= 0:
            return float("inf")
        radius_m = span_nm * 1e-9
        volume_m3 = (4.0 / 3.0) * math.pi * (radius_m ** 3)
        volume_liters = volume_m3 * 1000.0
        c_eff_m = 1.0 / (cls.AVOGADRO * volume_liters)
        return c_eff_m

    @classmethod
    def calculate_apparent_avidity(
        cls,
        kd_a: float,
        kd_b: float,
        linker_length_aa: int = 15,
        c_eff_override: Optional[float] = None,
    ) -> Dict[str, Any]:
        """
        Calculate apparent avidity KD_app and enhancement factor:
        KD_app = (KD_A * KD_B) / (C_eff + KD_B)
        Avidity Enhancement Factor beta = KD_A / KD_app
        """
        if kd_a <= 0 or kd_b <= 0:
            raise ValueError("KD values must be positive")

        if c_eff_override is not None:
            c_eff = c_eff_override
        else:
            lp = LinkerProperties("custom", linker_length_aa)
            span_nm = max(lp.root_mean_square_end_to_end_nm, 1.0)
            c_eff = cls.calculate_effective_concentration(span_nm)

        kd_app = (kd_a * kd_b) / (c_eff + kd_b)
        avidity_factor = kd_a / kd_app if kd_app > 0 else float("inf")
        cooperativity_index = math.log10(max(avidity_factor, 1.0))

        return {
            "monovalent_kd_a_m": kd_a,
            "monovalent_kd_b_m": kd_b,
            "effective_concentration_m": c_eff,
            "apparent_avidity_kd_m": kd_app,
            "avidity_enhancement_factor": avidity_factor,
            "cooperativity_index_log10": cooperativity_index,
            "is_avidity_driven": avidity_factor > 2.0,
        }


# ==============================================================================
# 3. T-CELL ENGAGEMENT & CYTOTOXICITY PREDICTOR
# ==============================================================================

class TCellEngagerPredictor:
    """Predicts immunological synapse formation, CD3 activation, and target cell lysis."""

    @staticmethod
    def predict_lysis(
        construct: BispecificConstruct,
        target_cell_density: float,  # receptors per tumor cell (e.g. 50,000 CD19)
        cd3_density: float = 40000.0,  # CD3 copies per T-cell (~40,000)
        effector_to_target_ratio: float = 5.0,  # E:T ratio (e.g. 5:1)
        antibody_conc_m: float = 1e-10,  # 100 pM
        incubation_time_hours: float = 24.0,
    ) -> Dict[str, Any]:
        """
        Simulate target cell lysis percentage based on synapse formation density.
        """
        kd_tumor = construct.arm_a.kd_molar
        kd_cd3 = construct.arm_b.kd_molar

        cell_conc_m = 1e6 * 1000 / 6.022e23
        tumor_receptor_m = target_cell_density * cell_conc_m
        tcell_cd3_m = cd3_density * cell_conc_m * effector_to_target_ratio

        ternary_res = TernaryComplexModel.solve_equilibrium(
            antibody_conc_m=antibody_conc_m,
            target_a_conc_m=tumor_receptor_m,
            target_b_conc_m=tcell_cd3_m,
            kd_a=kd_tumor,
            kd_b=kd_cd3,
            alpha=construct.cooperativity_alpha,
        )

        ternary_m = ternary_res["ternary_complex"]
        synapses_per_cell = (ternary_m / max(tumor_receptor_m, 1e-30)) * target_cell_density

        synapse_ec50 = 25.0
        hill_coeff = 2.2

        if synapses_per_cell <= 0:
            lysis_pct = 0.0
        else:
            time_factor = min(1.0, incubation_time_hours / 24.0)
            lysis_pct = (100.0 * (synapses_per_cell ** hill_coeff) /
                         ((synapse_ec50 ** hill_coeff) + (synapses_per_cell ** hill_coeff))) * time_factor

        crs_risk_index = min(10.0, (ternary_m / (kd_cd3 + 1e-30)) * 100.0)
        crs_category = "Low" if crs_risk_index < 2.5 else "Moderate" if crs_risk_index < 6.0 else "High"

        return {
            "construct_id": construct.construct_id,
            "antibody_conc_m": antibody_conc_m,
            "target_cell_density": target_cell_density,
            "cd3_density": cd3_density,
            "effector_to_target_ratio": effector_to_target_ratio,
            "synapses_per_tumor_cell": round(synapses_per_cell, 2),
            "specific_lysis_percentage": round(min(100.0, max(0.0, lysis_pct)), 2),
            "crs_risk_index": round(crs_risk_index, 2),
            "crs_risk_category": crs_category,
            "ternary_complex_molar": ternary_m,
        }


# ==============================================================================
# 4. SPR KINETICS & GLOBAL CURVE FITTER
# ==============================================================================

class SPRKineticsEngine:
    """SPR Sensorgram simulation and global parameter regression."""

    @staticmethod
    def simulate_sensorgram(
        ka: float,
        kd: float,
        rmax: float,
        concentration_m: float,
        time_points_s: List[float],
        association_duration_s: float,
        ri_bulk_refractive_index: float = 0.0,
    ) -> List[float]:
        """
        Simulate 1:1 Langmuir binding SPR sensorgram (RU vs time).
        """
        responses = []
        k_on_app = ka * concentration_m + kd
        req = (rmax * ka * concentration_m) / k_on_app if k_on_app > 0 else 0.0

        r_assoc_end = req * (1.0 - math.exp(-k_on_app * association_duration_s)) if k_on_app > 0 else 0.0

        for t in time_points_s:
            if t <= association_duration_s:
                if k_on_app > 0:
                    r_t = req * (1.0 - math.exp(-k_on_app * t)) + ri_bulk_refractive_index
                else:
                    r_t = 0.0
            else:
                t_diss = t - association_duration_s
                r_t = r_assoc_end * math.exp(-kd * t_diss)
            responses.append(r_t)

        return responses

    @classmethod
    def fit_global_1_1(
        cls,
        curves: List[SPRCurve],
        initial_ka: float = 1e5,
        initial_kd: float = 1e-3,
        initial_rmax: float = 100.0,
        iterations: int = 300,
        learning_rate: float = 0.05,
    ) -> Dict[str, Any]:
        """
        Global non-linear least-squares fitting of ka, kd, and Rmax across multiple concentrations
        using log-space parameterization for numerical stability.
        """
        if not curves:
            raise ValueError("No SPR curves provided for fitting")

        total_points = sum(len(c.time_s) for c in curves)
        if total_points == 0:
            raise ValueError("Curves contain zero time points")

        u = math.log(max(initial_ka, 1e2))
        v = math.log(max(initial_kd, 1e-6))
        w = math.log(max(initial_rmax, 1.0))

        best_rss = float("inf")
        best_params = (u, v, w)

        v_u, v_v, v_w = 0.0, 0.0, 0.0
        momentum = 0.8

        for it in range(iterations):
            ka = math.exp(u)
            kd = math.exp(v)
            rmax = math.exp(w)

            rss = 0.0
            grad_u = 0.0
            grad_v = 0.0
            grad_w = 0.0

            for curve in curves:
                conc = curve.concentration_m
                t_assoc = curve.association_duration_s
                k_app = ka * conc + kd
                req = (rmax * ka * conc) / k_app if k_app > 0 else 0.0
                r_end = req * (1.0 - math.exp(-k_app * t_assoc)) if k_app > 0 else 0.0

                for t, y_obs in zip(curve.time_s, curve.response_ru):
                    if t <= t_assoc:
                        exp_term = math.exp(-k_app * t)
                        y_pred = req * (1.0 - exp_term)
                        resid = y_pred - y_obs
                        rss += resid * resid

                        dreq_dka = (rmax * conc * kd) / (k_app ** 2) if k_app > 0 else 0.0
                        dreq_dkd = -(rmax * ka * conc) / (k_app ** 2) if k_app > 0 else 0.0
                        dreq_drmax = (ka * conc) / k_app if k_app > 0 else 0.0

                        dy_dka = dreq_dka * (1.0 - exp_term) + req * (conc * t * exp_term)
                        dy_dkd = dreq_dkd * (1.0 - exp_term) + req * (t * exp_term)
                        dy_drmax = dreq_drmax * (1.0 - exp_term)
                    else:
                        t_diss = t - t_assoc
                        exp_diss = math.exp(-kd * t_diss)
                        y_pred = r_end * exp_diss
                        resid = y_pred - y_obs
                        rss += resid * resid

                        dreq_dka = (rmax * conc * kd) / (k_app ** 2) if k_app > 0 else 0.0
                        dreq_dkd = -(rmax * ka * conc) / (k_app ** 2) if k_app > 0 else 0.0
                        dreq_drmax = (ka * conc) / k_app if k_app > 0 else 0.0

                        dr_end_dka = dreq_dka * (1.0 - math.exp(-k_app * t_assoc)) + req * (conc * t_assoc * math.exp(-k_app * t_assoc))
                        dr_end_dkd = dreq_dkd * (1.0 - math.exp(-k_app * t_assoc)) + req * (t_assoc * math.exp(-k_app * t_assoc))
                        dr_end_drmax = dreq_drmax * (1.0 - math.exp(-k_app * t_assoc))

                        dy_dka = dr_end_dka * exp_diss
                        dy_dkd = dr_end_dkd * exp_diss - r_end * t_diss * exp_diss
                        dy_drmax = dr_end_drmax * exp_diss

                    # Convert to log-derivatives: d/du = d/dka * ka
                    dy_du = dy_dka * ka
                    dy_dv = dy_dkd * kd
                    dy_dw = dy_drmax * rmax

                    grad_u += 2.0 * resid * dy_du
                    grad_v += 2.0 * resid * dy_dv
                    grad_w += 2.0 * resid * dy_dw

            if rss < best_rss:
                best_rss = rss
                best_params = (u, v, w)

            # Gradient clipping
            grad_norm = math.sqrt(grad_u ** 2 + grad_v ** 2 + grad_w ** 2) / max(total_points, 1)
            scale = min(1.0, 10.0 / max(grad_norm, 1e-12))

            lr = (learning_rate / (1.0 + 0.01 * it)) * scale
            v_u = momentum * v_u + lr * (grad_u / total_points)
            v_v = momentum * v_v + lr * (grad_v / total_points)
            v_w = momentum * v_w + lr * (grad_w / total_points)

            u = max(math.log(1e2), min(math.log(1e8), u - v_u))
            v = max(math.log(1e-7), min(math.log(1e1), v - v_v))
            w = max(math.log(1.0), min(math.log(1e4), w - v_w))

        u_opt, v_opt, w_opt = best_params
        ka_opt = math.exp(u_opt)
        kd_opt = math.exp(v_opt)
        rmax_opt = math.exp(w_opt)
        kd_molar = kd_opt / ka_opt if ka_opt > 0 else float("inf")

        all_obs = [y for c in curves for y in c.response_ru]
        mean_obs = sum(all_obs) / len(all_obs) if all_obs else 0.0
        sst = sum((y - mean_obs) ** 2 for y in all_obs)
        r_squared = max(0.0, 1.0 - (best_rss / sst)) if sst > 0 else 1.0
        chi_squared = best_rss / max(1, total_points - 3)

        return {
            "ka": ka_opt,
            "kd": kd_opt,
            "KD_molar": kd_molar,
            "Rmax_RU": rmax_opt,
            "total_points": total_points,
            "rss": best_rss,
            "chi_squared": chi_squared,
            "r_squared": round(r_squared, 4),
            "fitting_quality": "Excellent" if r_squared >= 0.98 else "Good" if r_squared >= 0.90 else "Suboptimal",
        }


# ==============================================================================
# 5. LINKER MECHANICS & STRUCTURAL OPTIMIZATION
# ==============================================================================

class LinkerOptimizer:
    """Evaluates linker length and biophysical strain for dual epitope bridging."""

    @classmethod
    def evaluate_linker(
        cls,
        linker_length_aa: int,
        target_epitope_distance_nm: float,
        kd_a: float,
        kd_b: float,
    ) -> Dict[str, Any]:
        """
        Evaluate structural feasibility and bridging efficiency of a linker.
        """
        lp = LinkerProperties("G4S_series", linker_length_aa)
        contour_len = lp.contour_length_nm
        rms_span = lp.root_mean_square_end_to_end_nm

        if target_epitope_distance_nm <= 0:
            target_epitope_distance_nm = 4.5

        is_too_short = contour_len < target_epitope_distance_nm
        extension_ratio = target_epitope_distance_nm / max(contour_len, 0.1)

        if extension_ratio >= 1.0:
            strain_penalty = 100.0
            geometric_feasibility = "Infeasible (Too Short)"
            bridge_score = 0.0
        elif extension_ratio > 0.85:
            strain_penalty = 25.0
            geometric_feasibility = "High Strain"
            bridge_score = 0.35
        elif extension_ratio < 0.2:
            strain_penalty = 5.0
            geometric_feasibility = "Excessive Flexibility / Low C_eff"
            bridge_score = 0.65
        else:
            strain_penalty = 0.0
            geometric_feasibility = "Optimal Geometry"
            bridge_score = 0.95

        c_eff = AvidityEngine.calculate_effective_concentration(max(rms_span, 1.0))
        apparent_avidity = AvidityEngine.calculate_apparent_avidity(kd_a, kd_b, linker_length_aa)

        composite_score = bridge_score * apparent_avidity["cooperativity_index_log10"]

        return {
            "linker_length_aa": linker_length_aa,
            "target_epitope_distance_nm": target_epitope_distance_nm,
            "contour_length_nm": round(contour_len, 2),
            "rms_end_to_end_nm": round(rms_span, 2),
            "extension_ratio": round(extension_ratio, 3),
            "geometric_feasibility": geometric_feasibility,
            "strain_penalty_kcal_mol": strain_penalty,
            "effective_concentration_m": c_eff,
            "apparent_avidity_kd_m": apparent_avidity["apparent_avidity_kd_m"],
            "composite_optimization_score": round(composite_score, 4),
        }

    @classmethod
    def find_optimal_linker(
        cls,
        target_epitope_distance_nm: float,
        kd_a: float,
        kd_b: float,
        candidate_lengths: Optional[List[int]] = None,
    ) -> Dict[str, Any]:
        """Rank candidate linker lengths to find the biophysical optimum."""
        if candidate_lengths is None:
            candidate_lengths = [5, 10, 15, 20, 25, 30, 35, 40]

        evaluations = [
            cls.evaluate_linker(n, target_epitope_distance_nm, kd_a, kd_b)
            for n in candidate_lengths
        ]

        evaluations.sort(key=lambda x: x["composite_optimization_score"], reverse=True)
        best = evaluations[0] if evaluations else None

        return {
            "target_epitope_distance_nm": target_epitope_distance_nm,
            "optimal_linker_length_aa": best["linker_length_aa"] if best else None,
            "best_evaluation": best,
            "ranked_candidates": evaluations,
        }
