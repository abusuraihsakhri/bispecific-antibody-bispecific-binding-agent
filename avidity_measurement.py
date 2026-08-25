"""Bispecific Avidity Measurement: avidity coefficient calculation and cooperative binding analysis."""
from typing import Dict, Any, List, Optional
from dataclasses import dataclass, field
import math


@dataclass
class BindingMeasurement:
    target_a_conc_m: float
    target_b_conc_m: float
    response_ru: float
    time_s: float


class AvidityMeasurement:
    """Calculate avidity coefficients from bispecific binding data."""

    def __init__(self):
        self._measurements: List[BindingMeasurement] = []

    def add_measurement(self, m: BindingMeasurement) -> None:
        self._measurements.append(m)

    def calculate_avidity_coefficient(
        self, kd_mono_a: float, kd_mono_b: float, measured_response: float,
        max_theoretical_response: float
    ) -> Dict[str, Any]:
        """Avidity coefficient = measured response / monovalent expectation."""
        mono_expectation = max_theoretical_response * 0.5
        avidity_coeff = measured_response / max(mono_expectation, 1e-10)
        return {
            "avidity_coefficient": round(avidity_coeff, 4),
            "measured_response_ru": round(measured_response, 2),
            "mono_expectation_ru": round(mono_expectation, 2),
            "enhancement_factor": round(avidity_coeff / 1.0, 4),
            "is_multivalent": avidity_coeff > 1.2,
        }

    def cooperative_binding_analysis(
        self, measurements: Optional[List[BindingMeasurement]] = None
    ) -> Dict[str, Any]:
        """Detect cooperative binding from concentration-response curves."""
        data = measurements or self._measurements
        if len(data) < 4:
            return {"status": "insufficient_data"}

        responses = [m.response_ru for m in data]
        concentrations = [m.target_a_conc_m for m in data]

        min_r, max_r = min(responses), max(responses)
        range_r = max_r - min_r if max_r > min_r else 1e-10

        normalized = [(r - min_r) / range_r for r in responses]

        hill_slopes = []
        for i in range(1, len(normalized) - 1):
            if concentrations[i] > 0 and concentrations[i - 1] > 0:
                if normalized[i] > 0 and normalized[i - 1] > 0:
                    log_ratio = math.log(normalized[i] / max(normalized[i - 1], 1e-10))
                    log_conc_ratio = math.log(concentrations[i] / max(concentrations[i - 1], 1e-10))
                    if abs(log_conc_ratio) > 1e-10:
                        hill_slopes.append(log_ratio / log_conc_ratio)

        avg_hill = sum(hill_slopes) / len(hill_slopes) if hill_slopes else 1.0

        if avg_hill > 1.5:
            cooperativity = "positive"
        elif avg_hill < 0.5:
            cooperativity = "negative"
        else:
            cooperativity = "non_cooperative"

        ec50_est = self._estimate_ec50(concentrations, normalized)

        return {
            "hill_coefficient": round(avg_hill, 4),
            "cooperativity": cooperativity,
            "ec50_estimate": ec50_est,
            "num_measurements": len(data),
        }

    def _estimate_ec50(self, concentrations: List[float], normalized: List[float]) -> float:
        for i in range(len(normalized) - 1):
            if normalized[i] <= 0.5 <= normalized[i + 1]:
                if normalized[i + 1] - normalized[i] > 1e-10:
                    frac = (0.5 - normalized[i]) / (normalized[i + 1] - normalized[i])
                    return concentrations[i] * (concentrations[i + 1] / max(concentrations[i], 1e-10)) ** frac
        return concentrations[len(concentrations) // 2] if concentrations else 0.0
