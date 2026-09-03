"""
Unit Test Suite for Bispecific Antibody Binding Agent & Kinetics Engine
========================================================================
Comprehensive verification across:
  - Ternary complex equilibrium & bell-shaped hook effect mechanics
  - Avidity enhancement and effective local concentration C_eff
  - T-cell immunological synapse formation and cytotoxicity prediction
  - Global SPR kinetics multi-curve fitting
  - Polymer physics linker mechanics optimization
  - Error handling and boundary conditions
"""

import math
import unittest
from bispecific_antibody import (
    AvidityEngine,
    BindingArm,
    BispecificConstruct,
    LinkerOptimizer,
    LinkerProperties,
    SPRCurve,
    SPRKineticsEngine,
    TCellEngagerPredictor,
    TernaryComplexModel,
)


class TestTernaryComplexModel(unittest.TestCase):
    """Test ternary complex equilibrium calculations and hook effect."""

    def test_zero_antibody_concentration(self):
        res = TernaryComplexModel.solve_equilibrium(
            antibody_conc_m=0.0,
            target_a_conc_m=1e-8,
            target_b_conc_m=1e-8,
            kd_a=1e-9,
            kd_b=1e-9,
        )
        self.assertEqual(res["ternary_complex"], 0.0)
        self.assertEqual(res["binary_ab_a"], 0.0)
        self.assertEqual(res["free_ab"], 0.0)

    def test_zero_target_concentration(self):
        res = TernaryComplexModel.solve_equilibrium(
            antibody_conc_m=1e-9,
            target_a_conc_m=0.0,
            target_b_conc_m=1e-8,
            kd_a=1e-9,
            kd_b=1e-9,
        )
        self.assertEqual(res["ternary_complex"], 0.0)
        self.assertEqual(res["binary_ab_a"], 0.0)

    def test_symmetric_binding_formation(self):
        res = TernaryComplexModel.solve_equilibrium(
            antibody_conc_m=1e-9,
            target_a_conc_m=1e-8,
            target_b_conc_m=1e-8,
            kd_a=1e-9,
            kd_b=1e-9,
            alpha=1.0,
        )
        self.assertGreater(res["ternary_complex"], 0.0)
        self.assertAlmostEqual(res["binary_ab_a"], res["binary_ab_b"], places=12)

    def test_asymmetric_affinities(self):
        res = TernaryComplexModel.solve_equilibrium(
            antibody_conc_m=1e-8,
            target_a_conc_m=1e-8,
            target_b_conc_m=1e-8,
            kd_a=1e-10,  # High affinity
            kd_b=1e-8,   # Moderate affinity
        )
        self.assertGreater(res["binary_ab_a"], res["binary_ab_b"])
        self.assertGreater(res["ternary_complex"], 0.0)

    def test_positive_cooperativity(self):
        res_neutral = TernaryComplexModel.solve_equilibrium(
            antibody_conc_m=1e-9,
            target_a_conc_m=1e-8,
            target_b_conc_m=1e-8,
            kd_a=1e-9,
            kd_b=1e-9,
            alpha=1.0,
        )
        res_synergy = TernaryComplexModel.solve_equilibrium(
            antibody_conc_m=1e-9,
            target_a_conc_m=1e-8,
            target_b_conc_m=1e-8,
            kd_a=1e-9,
            kd_b=1e-9,
            alpha=10.0,  # 10x positive cooperativity
        )
        self.assertGreater(res_synergy["ternary_complex"], res_neutral["ternary_complex"])

    def test_negative_cooperativity(self):
        res_neutral = TernaryComplexModel.solve_equilibrium(
            antibody_conc_m=1e-9,
            target_a_conc_m=1e-8,
            target_b_conc_m=1e-8,
            kd_a=1e-9,
            kd_b=1e-9,
            alpha=1.0,
        )
        res_antagonism = TernaryComplexModel.solve_equilibrium(
            antibody_conc_m=1e-9,
            target_a_conc_m=1e-8,
            target_b_conc_m=1e-8,
            kd_a=1e-9,
            kd_b=1e-9,
            alpha=0.1,  # Negative cooperativity
        )
        self.assertLess(res_antagonism["ternary_complex"], res_neutral["ternary_complex"])

    def test_hook_effect_titration(self):
        res = TernaryComplexModel.compute_dose_response_curve(
            target_a_conc_m=1e-8,
            target_b_conc_m=1e-8,
            kd_a=1e-9,
            kd_b=1e-9,
        )
        self.assertTrue(res["hook_effect_present"])
        self.assertGreater(res["max_ternary_complex_m"], 0.0)
        self.assertGreater(res["optimal_antibody_conc_m"], 0.0)

        # High antibody concentration should yield lower ternary complex than the peak
        highest_ab_pt = res["dose_response_points"][-1]
        self.assertLess(highest_ab_pt["ternary_complex_m"], res["max_ternary_complex_m"])

    def test_invalid_negative_concentrations(self):
        with self.assertRaises(ValueError):
            TernaryComplexModel.solve_equilibrium(
                antibody_conc_m=-1e-9,
                target_a_conc_m=1e-8,
                target_b_conc_m=1e-8,
                kd_a=1e-9,
                kd_b=1e-9,
            )

    def test_invalid_negative_kd(self):
        with self.assertRaises(ValueError):
            TernaryComplexModel.solve_equilibrium(
                antibody_conc_m=1e-9,
                target_a_conc_m=1e-8,
                target_b_conc_m=1e-8,
                kd_a=-1e-9,
                kd_b=1e-9,
            )


class TestAvidityEngine(unittest.TestCase):
    """Test avidity and effective concentration biophysics."""

    def test_effective_concentration_scaling(self):
        c_eff_small = AvidityEngine.calculate_effective_concentration(span_nm=3.0)
        c_eff_large = AvidityEngine.calculate_effective_concentration(span_nm=10.0)
        # Shorter span -> smaller search volume -> higher effective concentration
        self.assertGreater(c_eff_small, c_eff_large)

    def test_avidity_enhancement_calculation(self):
        res = AvidityEngine.calculate_apparent_avidity(
            kd_a=1e-8,
            kd_b=1e-8,
            linker_length_aa=15,
        )
        self.assertLess(res["apparent_avidity_kd_m"], res["monovalent_kd_a_m"])
        self.assertGreater(res["avidity_enhancement_factor"], 1.0)
        self.assertTrue(res["is_avidity_driven"])

    def test_avidity_with_override(self):
        c_eff_explicit = 1e-4  # 100 uM
        res = AvidityEngine.calculate_apparent_avidity(
            kd_a=1e-6,
            kd_b=1e-6,
            c_eff_override=c_eff_explicit,
        )
        self.assertEqual(res["effective_concentration_m"], c_eff_explicit)
        expected_kd = (1e-6 * 1e-6) / (1e-4 + 1e-6)
        self.assertAlmostEqual(res["apparent_avidity_kd_m"], expected_kd, places=12)

    def test_zero_or_negative_kd_error(self):
        with self.assertRaises(ValueError):
            AvidityEngine.calculate_apparent_avidity(kd_a=0.0, kd_b=1e-8)


class TestLinkerMechanics(unittest.TestCase):
    """Test polymer physics and linker mechanics."""

    def test_linker_contour_length(self):
        lp = LinkerProperties(sequence_or_type="(G4S)3", length_aa=15)
        self.assertAlmostEqual(lp.contour_length_nm, 15 * 0.38, places=4)

    def test_wlc_rms_end_to_end(self):
        lp = LinkerProperties(sequence_or_type="(G4S)3", length_aa=15)
        rms = lp.root_mean_square_end_to_end_nm
        self.assertGreater(rms, 0.0)
        self.assertLess(rms, lp.contour_length_nm)

    def test_linker_too_short_evaluation(self):
        res = LinkerOptimizer.evaluate_linker(
            linker_length_aa=5,  # Contour ~ 1.9 nm
            target_epitope_distance_nm=5.0,  # Requires >= 5 nm
            kd_a=1e-9,
            kd_b=1e-8,
        )
        self.assertEqual(res["geometric_feasibility"], "Infeasible (Too Short)")
        self.assertEqual(res["composite_optimization_score"], 0.0)

    def test_find_optimal_linker_ranking(self):
        res = LinkerOptimizer.find_optimal_linker(
            target_epitope_distance_nm=4.5,
            kd_a=1e-9,
            kd_b=1e-8,
            candidate_lengths=[5, 15, 30, 50],
        )
        self.assertIsNotNone(res["optimal_linker_length_aa"])
        self.assertGreater(res["optimal_linker_length_aa"], 5)
        self.assertLessEqual(res["optimal_linker_length_aa"], 30)


class TestTCellEngagerPredictor(unittest.TestCase):
    """Test immunological synapse and target cell cytotoxicity."""

    def setUp(self):
        self.construct = BispecificConstruct(
            construct_id="BiTE-CD19xCD3",
            name="Blinatumomab-like",
            arm_a=BindingArm("CD19", ka=1e5, kd=1e-4),
            arm_b=BindingArm("CD3", ka=1e5, kd=1e-3),
            linker=LinkerProperties("(G4S)3", length_aa=15),
        )

    def test_zero_antibody_lysis(self):
        res = TCellEngagerPredictor.predict_lysis(
            construct=self.construct,
            target_cell_density=50000.0,
            antibody_conc_m=0.0,
        )
        self.assertEqual(res["specific_lysis_percentage"], 0.0)
        self.assertEqual(res["synapses_per_tumor_cell"], 0.0)

    def test_potent_lysis_at_nanomolar(self):
        res = TCellEngagerPredictor.predict_lysis(
            construct=self.construct,
            target_cell_density=50000.0,
            antibody_conc_m=1e-9,  # 1 nM
            effector_to_target_ratio=10.0,
            incubation_time_hours=24.0,
        )
        self.assertGreater(res["specific_lysis_percentage"], 50.0)
        self.assertGreater(res["synapses_per_tumor_cell"], 10.0)

    def test_crs_risk_levels(self):
        res_low = TCellEngagerPredictor.predict_lysis(
            construct=self.construct,
            target_cell_density=50000.0,
            antibody_conc_m=1e-12,
        )
        res_high = TCellEngagerPredictor.predict_lysis(
            construct=self.construct,
            target_cell_density=50000.0,
            antibody_conc_m=1e-6,
        )
        self.assertLessEqual(res_low["crs_risk_index"], res_high["crs_risk_index"])


class TestSPRKineticsEngine(unittest.TestCase):
    """Test SPR sensorgram simulation and non-linear regression."""

    def test_sensorgram_shape(self):
        ka, kd, rmax = 1e5, 1e-3, 100.0
        t_points = [0.0, 30.0, 60.0, 90.0, 120.0]
        assoc_t = 60.0
        responses = SPRKineticsEngine.simulate_sensorgram(
            ka=ka, kd=kd, rmax=rmax,
            concentration_m=1e-8,
            time_points_s=t_points,
            association_duration_s=assoc_t,
        )
        self.assertEqual(len(responses), len(t_points))
        # Association: R(30) > R(0)
        self.assertGreater(responses[1], responses[0])
        # Dissociation: R(120) < R(60)
        self.assertLess(responses[4], responses[2])

    def test_global_curve_fitting(self):
        true_ka, true_kd, true_rmax = 2e5, 2e-3, 100.0
        concs = [1e-8, 3e-8, 1e-7]
        t_points = [float(i) for i in range(0, 120, 5)]

        curves = []
        for c in concs:
            resp = SPRKineticsEngine.simulate_sensorgram(true_ka, true_kd, true_rmax, c, t_points, 60.0)
            curves.append(SPRCurve(c, t_points, resp, 60.0))

        fit = SPRKineticsEngine.fit_global_1_1(curves, initial_ka=1e5, initial_kd=1e-3, initial_rmax=80.0)
        self.assertGreaterEqual(fit["r_squared"], 0.95)
        self.assertGreater(fit["ka"], 0)
        self.assertGreater(fit["kd"], 0)
        self.assertGreater(fit["KD_molar"], 0)

    def test_empty_curves_error(self):
        with self.assertRaises(ValueError):
            SPRKineticsEngine.fit_global_1_1([])

    def test_zero_timepoints_error(self):
        empty_curve = SPRCurve(1e-8, [], [], 60.0)
        with self.assertRaises(ValueError):
            SPRKineticsEngine.fit_global_1_1([empty_curve])


class TestCLIExecution(unittest.TestCase):
    """Test CLI commands execution."""

    def test_cli_ternary_json(self):
        import io
        import json
        from contextlib import redirect_stdout
        from cli import main

        buf = io.StringIO()
        with redirect_stdout(buf):
            ret = main(["ternary", "--ab-conc", "1e-9", "--json"])
        self.assertEqual(ret, 0)
        data = json.loads(buf.getvalue())
        self.assertIn("ternary_complex", data)

    def test_cli_titration_json(self):
        import io
        import json
        from contextlib import redirect_stdout
        from cli import main

        buf = io.StringIO()
        with redirect_stdout(buf):
            ret = main(["titration", "--target-a-conc", "1e-8", "--json"])
        self.assertEqual(ret, 0)
        data = json.loads(buf.getvalue())
        self.assertTrue(data["hook_effect_present"])

    def test_cli_avidity_json(self):
        import io
        import json
        from contextlib import redirect_stdout
        from cli import main

        buf = io.StringIO()
        with redirect_stdout(buf):
            ret = main(["avidity", "--kd-a", "1e-8", "--kd-b", "1e-8", "--json"])
        self.assertEqual(ret, 0)
        data = json.loads(buf.getvalue())
        self.assertIn("apparent_avidity_kd_m", data)

    def test_cli_tcell_lysis_json(self):
        import io
        import json
        from contextlib import redirect_stdout
        from cli import main

        buf = io.StringIO()
        with redirect_stdout(buf):
            ret = main(["tcell-lysis", "--tumor-density", "50000", "--json"])
        self.assertEqual(ret, 0)
        data = json.loads(buf.getvalue())
        self.assertIn("specific_lysis_percentage", data)

    def test_cli_optimize_linker_json(self):
        import io
        import json
        from contextlib import redirect_stdout
        from cli import main

        buf = io.StringIO()
        with redirect_stdout(buf):
            ret = main(["optimize-linker", "--distance", "4.5", "--json"])
        self.assertEqual(ret, 0)
        data = json.loads(buf.getvalue())
        self.assertIn("optimal_linker_length_aa", data)

    def test_cli_fit_spr_json(self):
        import io
        import json
        from contextlib import redirect_stdout
        from cli import main

        buf = io.StringIO()
        with redirect_stdout(buf):
            ret = main(["fit-spr", "--json"])
        self.assertEqual(ret, 0)
        data = json.loads(buf.getvalue())
        self.assertIn("r_squared", data)

    def test_cli_batch_csv_and_json(self):
        import io
        import json
        import os
        import tempfile
        from contextlib import redirect_stdout
        from cli import main

        sample_csv_content = (
            "construct_id,name,target_a,target_b,ab_conc_m,target_a_conc_m,target_b_conc_m,kd_a_m,kd_b_m,alpha,linker_length_aa,tumor_density,cd3_density,et_ratio\n"
            "BSAB-TEST,Test-BiTE,CD19,CD3,1.0e-10,5.0e-9,5.0e-9,1.5e-9,1.0e-7,1.0,15,50000,40000,5.0\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False, encoding="utf-8") as f_in:
            f_in.write(sample_csv_content)
            in_path = f_in.name

        out_path = in_path.replace(".csv", "_out.csv")

        try:
            # Test CSV output
            ret_csv = main(["batch", "-i", in_path, "-o", out_path])
            self.assertEqual(ret_csv, 0)
            self.assertTrue(os.path.exists(out_path))

            # Test JSON output
            buf = io.StringIO()
            with redirect_stdout(buf):
                ret_json = main(["batch", "-i", in_path, "--json"])
            self.assertEqual(ret_json, 0)
            data = json.loads(buf.getvalue())
            self.assertIsInstance(data, list)
            self.assertEqual(len(data), 1)
            self.assertEqual(data[0]["construct_id"], "BSAB-TEST")
            self.assertIn("ternary_complex_m", data[0])
            self.assertIn("apparent_avidity_kd_m", data[0])
        finally:
            if os.path.exists(in_path):
                os.remove(in_path)
            if os.path.exists(out_path):
                os.remove(out_path)


if __name__ == "__main__":
    unittest.main()

