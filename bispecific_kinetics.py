"""Bispecific Antibody Binding Kinetics: avidity modeling, cross-linking analysis, dual-target binding."""
from typing import Dict, Any, List, Optional
from dataclasses import dataclass, field
import math


@dataclass
class BispecificDesign:
    antibody_id: str
    target_a_name: str
    target_b_name: str
    ka_a: float
    kd_a: float
    ka_b: float
    kd_b: float
    linker_length_aa: int = 15
    valency_a: int = 1
    valency_b: int = 1


class BispecificBindingKinetics:
    """Modeling bispecific antibody binding to two targets simultaneously."""

    def calculate_monovalent_kd(self, ka: float, kd: float) -> float:
        """Simple K_D = k_off / k_on."""
        return kd / ka if ka > 0 else float('inf')

    def calculate_avidity_kd(
        self, kd_a: float, kd_b: float, effective_concentration: float = 1e-9
    ) -> Dict[str, Any]:
        """Effective avidity K_D from both binding arms."""
        kd_mono = math.sqrt(kd_a * kd_b)
        if effective_concentration > 0:
            avidity_enhancement = (kd_a + kd_b) / (2 * effective_concentration)
        else:
            avidity_enhancement = 1.0
        effective_kd = kd_mono / max(avidity_enhancement, 1.0)
        return {
            "mono_kd_a": kd_a, "mono_kd_b": kd_b,
            "geometric_mean_kd": kd_mono,
            "avidity_enhancement": round(avidity_enhancement, 4),
            "effective_avidity_kd": effective_kd,
        }

    def cross_linking_efficiency(
        self, design: BispecificDesign, target_a_conc: float, target_b_conc: float
    ) -> Dict[str, Any]:
        """Estimate efficiency of simultaneous target cross-linking."""
        ka_eff = math.sqrt(design.ka_a * design.ka_b)
        kd_eff = math.sqrt(design.kd_a * design.kd_b)
        avail_a = min(target_a_conc / (target_a_conc + kd_eff), 1.0)
        avail_b = min(target_b_conc / (target_b_conc + kd_eff), 1.0)
        cross_link_prob = avail_a * avail_b

        if design.linker_length_aa < 10:
            steric_risk = "high"
        elif design.linker_length_aa < 20:
            steric_risk = "moderate"
        else:
            steric_risk = "low"

        return {
            "ka_effective": round(ka_eff, 4),
            "kd_effective": round(kd_eff, 6),
            "target_a_occupancy": round(avail_a, 4),
            "target_b_occupancy": round(avail_b, 4),
            "cross_link_probability": round(cross_link_prob, 4),
            "steric_risk": steric_risk,
            "linker_length_aa": design.linker_length_aa,
        }

    def sequential_binding_kinetics(
        self, design: BispecificDesign, target_a_conc: float, target_b_conc: float
    ) -> Dict[str, Any]:
        """Model sequential binding: first arm binds, then second arm."""
        kd_a = self.calculate_monovalent_kd(design.ka_a, design.kd_a)
        kd_b = self.calculate_monovalent_kd(design.ka_b, design.kd_b)

        occupancy_a_first = target_a_conc / (target_a_conc + kd_a)
        occupancy_b_second = target_b_conc / (target_b_conc + kd_b)
        sequential_prob = occupancy_a_first * occupancy_b_second

        occupancy_b_first = target_b_conc / (target_b_conc + kd_b)
        occupancy_a_second = target_a_conc / (target_a_conc + kd_a)
        reverse_seq_prob = occupancy_b_first * occupancy_a_second

        return {
            "pathway_1_a_then_b": round(sequential_prob, 4),
            "pathway_2_b_then_a": round(reverse_seq_prob, 4),
            "dominant_pathway": "A_then_B" if sequential_prob >= reverse_seq_prob else "B_then_A",
            "overall_bivalent_occupancy": round(max(sequential_prob, reverse_seq_prob), 4),
        }

    def design_comparison(
        self, designs: List[BispecificDesign], target_a_conc: float = 10e-9, target_b_conc: float = 10e-9
    ) -> Dict[str, Any]:
        """Compare multiple bispecific designs."""
        results = []
        for d in designs:
            avidity = self.calculate_avidity_kd(d.kd_a, d.kd_b)
            cross_link = self.cross_link_efficiency(d, target_a_conc, target_b_conc)
            results.append({
                "antibody_id": d.antibody_id,
                "targets": f"{d.target_a_name} x {d.target_b_name}",
                "effective_avidity_kd": avidity["effective_avidity_kd"],
                "cross_link_probability": cross_link["cross_link_probability"],
                "steric_risk": cross_link["steric_risk"],
            })
        results.sort(key=lambda x: x["effective_avidity_kd"])
        for i, r in enumerate(results):
            r["rank"] = i + 1
        return {"ranked_designs": results, "best": results[0] if results else None}
