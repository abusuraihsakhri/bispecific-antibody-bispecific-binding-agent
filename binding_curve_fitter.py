"""Binding Curve Fitting: global fitting of multi-state binding models with parameter linkage."""
from typing import Dict, Any, List, Optional
from dataclasses import dataclass, field
import math


@dataclass
class BindingCurve:
    time_s: List[float]
    response_ru: List[float]
    concentration_m: float
    buffer_name: str = "default"


class GlobalBindingFitter:
    """Global fitting engine for multi-curve SPR binding data."""

    def __init__(self):
        self._curves: List[BindingCurve] = []

    def add_curve(self, curve: BindingCurve) -> None:
        self._curves.append(curve)

    def _model_1_1(self, t: float, ka: float, kd: float, rmax: float, c: float) -> float:
        k = ka * c + kd
        if k == 0:
            return 0
        req = rmax * c / (c + kd / ka) if ka > 0 else 0
        return req * (1 - math.exp(-k * t))

    def _model_bivalent(self, t: float, ka1: float, kd1: float, ka2: float, kd2: float,
                        rmax: float, c: float, f1: float = 0.5) -> float:
        r1 = self._model_1_1(t, ka1, kd1, rmax * f1, c)
        r2 = self._model_1_1(t, ka2, kd2, rmax * (1 - f1), c)
        return r1 + r2

    def fit_global_1_1(self, linkage: Optional[Dict[str, str]] = None) -> Dict[str, Any]:
        """Global fit all curves with shared or per-curve parameters."""
        if not self._curves:
            return {"status": "no_curves"}

        ka_est, kd_est, rmax_est = 1e5, 1e-3, 200.0
        lr = 1e-7

        for _ in range(100):
            total_rss = 0
            grad_ka = 0.0
            grad_kd = 0.0
            grad_rmax = 0.0

            for curve in self._curves:
                for t, r_obs in zip(curve.time_s, curve.response_ru):
                    pred = self._model_1_1(t, ka_est, kd_est, rmax_est, curve.concentration_m)
                    residual = r_obs - pred
                    total_rss += residual ** 2

                    k = ka_est * curve.concentration_m + kd_est
                    if k > 0 and ka_est > 0:
                        req = rmax_est * curve.concentration_m / (curve.concentration_m + kd_est / ka_est)
                        dp_dka = req * (1 - math.exp(-k * t)) * curve.concentration_m / ka_est
                        dp_dkd = -req * t * math.exp(-k * t)
                        dp_drmax = (curve.concentration_m / (curve.concentration_m + kd_est / ka_est)) * (1 - math.exp(-k * t))
                        grad_ka += -2 * residual * dp_dka
                        grad_kd += -2 * residual * dp_dkd
                        grad_rmax += -2 * residual * dp_drmax

            ka_est = max(ka_est - lr * grad_ka, 1e2)
            kd_est = max(kd_est - lr * grad_kd, 1e-6)
            rmax_est = max(rmax_est - lr * grad_rmax, 1.0)

        n_points = sum(len(c.time_s) for c in self._curves)
        rss = total_rss / max(n_points, 1)
        r_squared = max(0, 1 - rss / max(sum(c.response_ru[i] ** 2 for i in range(len(c.response_ru))) for c in self._curves), 0)

        return {
            "ka": round(ka_est, 4), "kd": round(kd_est, 6),
            "KD": round(kd_est / ka_est, 8) if ka_est > 0 else float('inf'),
            "rmax": round(rmax_est, 2),
            "rss": round(rss, 4), "r_squared": round(r_squared, 4),
            "num_curves": len(self._curves),
        }

    def residual_analysis(self, params: Dict[str, float]) -> Dict[str, Any]:
        """Analyze residuals across all curves for systematic deviations."""
        ka = params.get("ka", 1e5)
        kd = params.get("kd", 1e-3)
        rmax = params.get("rmax", 200.0)

        all_residuals = []
        curve_stats = []

        for curve in self._curves:
            residuals = []
            for t, r_obs in zip(curve.time_s, curve.response_ru):
                pred = self._model_1_1(t, ka, kd, rmax, curve.concentration_m)
                residuals.append(r_obs - pred)
            all_residuals.extend(residuals)

            mean_r = sum(residuals) / len(residuals) if residuals else 0
            max_abs_r = max(abs(r) for r in residuals) if residuals else 0
            curve_stats.append({
                "buffer": curve.buffer_name,
                "mean_residual": round(mean_r, 4),
                "max_abs_residual": round(max_abs_r, 4),
                "n_points": len(residuals),
            })

        overall_mean = sum(all_residuals) / len(all_residuals) if all_residuals else 0
        overall_std = math.sqrt(sum((r - overall_mean) ** 2 for r in all_residuals) / max(len(all_residuals) - 1, 1))

        return {
            "overall_mean_residual": round(overall_mean, 4),
            "overall_std_residual": round(overall_std, 4),
            "curve_details": curve_stats,
            "systematic_bias": abs(overall_mean) > 2 * overall_std,
        }
