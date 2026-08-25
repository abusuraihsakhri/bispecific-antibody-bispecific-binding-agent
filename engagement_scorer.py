"""Bivalent Target Engagement: avidity scoring, linker optimization, cross-linking probability."""
from typing import Dict, Any, List, Optional
from dataclasses import dataclass, field
import math


@dataclass
class LinkerDesign:
    length_aa: int
    composition: str
    flexibility: float
    hydrophobicity: float


class BivalentEngagementScorer:
    """Score bivalent target engagement efficiency and optimize linker design."""

    def score_engagement(
        self,
        kd_a: float, kd_b: float,
        target_a_conc: float, target_b_conc: float,
        linker: LinkerDesign,
    ) -> Dict[str, Any]:
        """Composite engagement score from binding affinities and geometry."""
        mono_score_a = min(1.0, 1e-6 / max(kd_a, 1e-15))
        mono_score_b = min(1.0, 1e-6 / max(kd_b, 1e-15))
        binding_score = math.sqrt(mono_score_a * mono_score_b)

        avail_a = target_a_conc / (target_a_conc + kd_a) if (target_a_conc + kd_a) > 0 else 0
        avail_b = target_b_conc / (target_b_conc + kd_b) if (target_b_conc + kd_b) > 0 else 0
        occupancy_score = avail_a * avail_b

        if linker.length_aa < 5:
            geometry_score = 0.2
        elif linker.length_aa < 12:
            geometry_score = 0.5
        elif linker.length_aa < 25:
            geometry_score = 0.9
        else:
            geometry_score = 0.7

        flexibility_bonus = 0.1 * (1 - linker.flexibility)
        geometry_score = min(1.0, geometry_score + flexibility_bonus)

        composite = 0.4 * binding_score + 0.3 * occupancy_score + 0.3 * geometry_score

        return {
            "composite_score": round(composite, 4),
            "binding_score": round(binding_score, 4),
            "occupancy_score": round(occupancy_score, 4),
            "geometry_score": round(geometry_score, 4),
            "cross_link_probability": round(occupancy_score * geometry_score, 4),
        }

    def optimize_linker(
        self, kd_a: float, kd_b: float, target_a_conc: float, target_b_conc: float,
        candidate_linkers: List[LinkerDesign],
    ) -> Dict[str, Any]:
        """Find optimal linker from candidates."""
        scored = []
        for linker in candidate_linkers:
            result = self.score_engagement(kd_a, kd_b, target_a_conc, target_b_conc, linker)
            scored.append({
                "length_aa": linker.length_aa,
                "composition": linker.composition,
                "score": result["composite_score"],
                "cross_link_probability": result["cross_link_probability"],
            })
        scored.sort(key=lambda x: x["score"], reverse=True)
        return {"ranked_linkers": scored, "best": scored[0] if scored else None}
