#!/usr/bin/env python3
"""
Command-Line Interface for Bispecific Antibody Binding Agent & Kinetics Engine
==============================================================================
Provides subcommands and interactive mode for:
  - Ternary complex equilibrium & bell-shaped hook effect
  - Avidity enhancement and effective local concentration C_eff
  - T-cell immunological synapse formation and cytotoxicity prediction
  - Global SPR kinetics multi-curve fitting
  - Structural linker mechanics optimization
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from typing import Any, Dict, List, Optional

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


def format_scientific(val: float) -> str:
    """Format floating point numbers in scientific notation or standard decimals."""
    if val == 0:
        return "0.00"
    if abs(val) < 1e-3 or abs(val) >= 1e5:
        return f"{val:.4e}"
    return f"{val:.4f}"


def cmd_ternary(args: argparse.Namespace) -> int:
    res = TernaryComplexModel.solve_equilibrium(
        antibody_conc_m=args.ab_conc,
        target_a_conc_m=args.target_a_conc,
        target_b_conc_m=args.target_b_conc,
        kd_a=args.kd_a,
        kd_b=args.kd_b,
        alpha=args.alpha,
    )
    if args.json:
        print(json.dumps(res, indent=2))
    else:
        print("=" * 60)
        print("  BISPECIFIC TERNARY COMPLEX EQUILIBRIUM ANALYSIS")
        print("=" * 60)
        print(f"Antibody Conc [Ab]0   : {format_scientific(args.ab_conc)} M")
        print(f"Target A Conc [R_A]0  : {format_scientific(args.target_a_conc)} M (KD = {format_scientific(args.kd_a)} M)")
        print(f"Target B Conc [R_B]0  : {format_scientific(args.target_b_conc)} M (KD = {format_scientific(args.kd_b)} M)")
        print(f"Cooperativity (Alpha) : {args.alpha:.2f}")
        print("-" * 60)
        print(f"Ternary Complex [A:Ab:B]: {format_scientific(res['ternary_complex'])} M")
        print(f"Target A Occupancy      : {res['ternary_fraction_of_target_a'] * 100:.2f}%")
        print(f"Target B Occupancy      : {res['ternary_fraction_of_target_b'] * 100:.2f}%")
        print(f"Binary Ab:TargetA       : {format_scientific(res['binary_ab_a'])} M")
        print(f"Binary Ab:TargetB       : {format_scientific(res['binary_ab_b'])} M")
        print(f"Free Unbound Antibody   : {format_scientific(res['free_ab'])} M")
        print("=" * 60)
    return 0


def cmd_titration(args: argparse.Namespace) -> int:
    res = TernaryComplexModel.compute_dose_response_curve(
        target_a_conc_m=args.target_a_conc,
        target_b_conc_m=args.target_b_conc,
        kd_a=args.kd_a,
        kd_b=args.kd_b,
        alpha=args.alpha,
    )
    if args.json:
        print(json.dumps(res, indent=2))
    else:
        print("=" * 65)
        print("  DOSE-RESPONSE & HOOK EFFECT TITRATION PROFILE")
        print("=" * 65)
        print(f"Peak Ternary Complex    : {format_scientific(res['max_ternary_complex_m'])} M")
        print(f"Optimal Antibody Conc   : {format_scientific(res['optimal_antibody_conc_m'])} M")
        print(f"Theoretical Optimum     : {format_scientific(res['theoretical_optimal_ab_m'])} M")
        print(f"Hook Effect Detected    : {'YES (Bell-shaped attenuation)' if res['hook_effect_present'] else 'NO'}")
        print("-" * 65)
        print(f"{'[Ab] (M)':<16} {'[Ternary] (M)':<16} {'[Ab:A] (M)':<16} {'[Ab:B] (M)':<16}")
        print("-" * 65)
        for pt in res["dose_response_points"][::6]:  # sample representative points
            print(f"{format_scientific(pt['antibody_conc_m']):<16} "
                  f"{format_scientific(pt['ternary_complex_m']):<16} "
                  f"{format_scientific(pt['binary_ab_a_m']):<16} "
                  f"{format_scientific(pt['binary_ab_b_m']):<16}")
        print("=" * 65)
    return 0


def cmd_avidity(args: argparse.Namespace) -> int:
    res = AvidityEngine.calculate_apparent_avidity(
        kd_a=args.kd_a,
        kd_b=args.kd_b,
        linker_length_aa=args.linker_len,
        c_eff_override=args.c_eff,
    )
    if args.json:
        print(json.dumps(res, indent=2))
    else:
        print("=" * 60)
        print("  BIVALENT AVIDITY & EFFECTIVE LOCAL CONCENTRATION")
        print("=" * 60)
        print(f"Arm A Monovalent KD     : {format_scientific(res['monovalent_kd_a_m'])} M")
        print(f"Arm B Monovalent KD     : {format_scientific(res['monovalent_kd_b_m'])} M")
        print(f"Effective Conc (C_eff)  : {format_scientific(res['effective_concentration_m'])} M")
        print(f"Apparent Avidity KD,app : {format_scientific(res['apparent_avidity_kd_m'])} M")
        print(f"Avidity Enhancement (β) : {res['avidity_enhancement_factor']:.2f}x")
        print(f"Avidity-Driven Binding  : {'YES' if res['is_avidity_driven'] else 'NO'}")
        print("=" * 60)
    return 0


def cmd_tcell_lysis(args: argparse.Namespace) -> int:
    construct = BispecificConstruct(
        construct_id=args.construct_id,
        name=args.name,
        arm_a=BindingArm(target_name=args.tumor_antigen, ka=1e5, kd=args.tumor_kd * 1e5),
        arm_b=BindingArm(target_name="CD3", ka=1e5, kd=args.cd3_kd * 1e5),
        linker=LinkerProperties(sequence_or_type="(G4S)3", length_aa=15),
        cooperativity_alpha=args.alpha,
    )
    res = TCellEngagerPredictor.predict_lysis(
        construct=construct,
        target_cell_density=args.tumor_density,
        cd3_density=args.cd3_density,
        effector_to_target_ratio=args.et_ratio,
        antibody_conc_m=args.ab_conc,
        incubation_time_hours=args.hours,
    )
    if args.json:
        print(json.dumps(res, indent=2))
    else:
        print("=" * 65)
        print("  T-CELL ENGAGER (BiTE) CYTOTOXICITY SIMULATION")
        print("=" * 65)
        print(f"Construct ID            : {res['construct_id']}")
        print(f"Tumor Antigen / CD3 KD  : {format_scientific(args.tumor_kd)} M / {format_scientific(args.cd3_kd)} M")
        print(f"Tumor Receptor Density  : {args.tumor_density:,.0f} receptors/cell")
        print(f"E:T Ratio / Incubation  : {args.et_ratio}:1 ({args.hours:.1f} hours)")
        print(f"Antibody Concentration  : {format_scientific(args.ab_conc)} M")
        print("-" * 65)
        print(f"Synapses per Tumor Cell : {res['synapses_per_tumor_cell']}")
        print(f"Specific Lysis Rate     : {res['specific_lysis_percentage']:.2f}%")
        print(f"CRS Risk Category       : {res['crs_risk_category']} (Index: {res['crs_risk_index']})")
        print("=" * 65)
    return 0


def cmd_fit_spr(args: argparse.Namespace) -> int:
    # Generate synthetic multi-cycle SPR curves for demonstration/testing
    t_points = [float(i) for i in range(0, 180, 2)]
    true_ka, true_kd, true_rmax = 2.5e5, 1.2e-3, 120.0
    concs = [5e-9, 15e-9, 50e-9, 150e-9]

    curves = []
    for c in concs:
        resp = SPRKineticsEngine.simulate_sensorgram(true_ka, true_kd, true_rmax, c, t_points, association_duration_s=60.0)
        curves.append(SPRCurve(concentration_m=c, time_s=t_points, response_ru=resp, association_duration_s=60.0))

    fit_res = SPRKineticsEngine.fit_global_1_1(curves, initial_ka=1e5, initial_kd=1e-3, initial_rmax=100.0)
    if args.json:
        print(json.dumps(fit_res, indent=2))
    else:
        print("=" * 60)
        print("  SPR GLOBAL 1:1 LANGMUIR KINETICS REGRESSION")
        print("=" * 60)
        print(f"Number of Curves Fitted : {len(curves)} concentrations")
        print(f"Total Datapoints        : {fit_res['total_points']}")
        print("-" * 60)
        print(f"Association Rate ka     : {format_scientific(fit_res['ka'])} M^-1 s^-1")
        print(f"Dissociation Rate kd    : {format_scientific(fit_res['kd'])} s^-1")
        print(f"Equilibrium KD          : {format_scientific(fit_res['KD_molar'])} M")
        print(f"Theoretical Rmax        : {fit_res['Rmax_RU']:.2f} RU")
        print(f"Goodness of Fit (R^2)   : {fit_res['r_squared']:.4f} ({fit_res['fitting_quality']})")
        print(f"Chi-Squared (χ^2)       : {fit_res['chi_squared']:.4f}")
        print("=" * 60)
    return 0


def cmd_optimize_linker(args: argparse.Namespace) -> int:
    res = LinkerOptimizer.find_optimal_linker(
        target_epitope_distance_nm=args.distance,
        kd_a=args.kd_a,
        kd_b=args.kd_b,
    )
    if args.json:
        print(json.dumps(res, indent=2))
    else:
        print("=" * 70)
        print("  POLYMER PHYSICS LINKER MECHANICS OPTIMIZATION")
        print("=" * 70)
        print(f"Target Epitope Distance : {args.distance:.2f} nm")
        print(f"Optimal Linker Length   : {res['optimal_linker_length_aa']} amino acids")
        print("-" * 70)
        print(f"{'Length (aa)':<12} {'Contour (nm)':<14} {'RMS Span (nm)':<15} {'Feasibility':<20} {'Score':<10}")
        print("-" * 70)
        for cand in res["ranked_candidates"]:
            print(f"{cand['linker_length_aa']:<12} "
                  f"{cand['contour_length_nm']:<14} "
                  f"{cand['rms_end_to_end_nm']:<15} "
                  f"{cand['geometric_feasibility']:<20} "
                  f"{cand['composite_optimization_score']:<10.4f}")
        print("=" * 70)
    return 0


def cmd_interactive() -> int:
    print("Bispecific Antibody Analytics Interactive REPL")
    print("Type 'help' for commands, 'exit' to quit.\n")
    while True:
        try:
            line = input("bsab> ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nExiting.")
            break
        if not line:
            continue
        if line.lower() in ("exit", "quit"):
            break
        if line.lower() == "help":
            print("Commands: ternary, titration, avidity, lysis, optimize, exit")
            continue

        parts = line.split()
        cmd = parts[0].lower()
        if cmd == "ternary":
            res = TernaryComplexModel.solve_equilibrium(1e-9, 10e-9, 10e-9, 1e-9, 1e-9)
            print(f"Ternary complex: {format_scientific(res['ternary_complex'])} M")
        elif cmd == "avidity":
            res = AvidityEngine.calculate_apparent_avidity(1e-9, 1e-9, 15)
            print(f"Apparent KD: {format_scientific(res['apparent_avidity_kd_m'])} M, enhancement: {res['avidity_enhancement_factor']:.1f}x")
        elif cmd == "lysis":
            c = BispecificConstruct("B01", "BiTE-CD19xCD3",
                                   BindingArm("CD19", 1e5, 1e-4),
                                   BindingArm("CD3", 1e5, 1e-3),
                                   LinkerProperties("(G4S)3", 15))
            res = TCellEngagerPredictor.predict_lysis(c, 50000.0, 40000.0, 5.0, 1e-10)
            print(f"Lysis: {res['specific_lysis_percentage']:.1f}%, CRS risk: {res['crs_risk_category']}")
        elif cmd == "optimize":
            res = LinkerOptimizer.find_optimal_linker(4.5, 1e-9, 1e-9)
            print(f"Optimal linker length: {res['optimal_linker_length_aa']} aa")
        else:
            print(f"Unknown command: {cmd}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="bispecific_antibody_cli",
        description="Bispecific Antibody Biophysics & Kinetics Analytics Platform",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: ternary
    p_tern = subparsers.add_parser("ternary", help="Calculate ternary complex equilibrium")
    p_tern.add_argument("--ab-conc", type=float, default=1e-9, help="Total antibody concentration (M)")
    p_tern.add_argument("--target-a-conc", type=float, default=1e-8, help="Target A concentration (M)")
    p_tern.add_argument("--target-b-conc", type=float, default=1e-8, help="Target B concentration (M)")
    p_tern.add_argument("--kd-a", type=float, default=1e-9, help="KD for Target A (M)")
    p_tern.add_argument("--kd-b", type=float, default=1e-8, help="KD for Target B (M)")
    p_tern.add_argument("--alpha", type=float, default=1.0, help="Cooperativity factor")
    p_tern.add_argument("--json", action="store_true", help="Output JSON")

    # Subcommand: titration
    p_titr = subparsers.add_parser("titration", help="Compute hook effect titration profile")
    p_titr.add_argument("--target-a-conc", type=float, default=1e-8, help="Target A concentration (M)")
    p_titr.add_argument("--target-b-conc", type=float, default=1e-8, help="Target B concentration (M)")
    p_titr.add_argument("--kd-a", type=float, default=1e-9, help="KD for Target A (M)")
    p_titr.add_argument("--kd-b", type=float, default=1e-8, help="KD for Target B (M)")
    p_titr.add_argument("--alpha", type=float, default=1.0, help="Cooperativity factor")
    p_titr.add_argument("--json", action="store_true", help="Output JSON")

    # Subcommand: avidity
    p_avid = subparsers.add_parser("avidity", help="Calculate effective concentration & avidity")
    p_avid.add_argument("--kd-a", type=float, default=1e-8, help="KD for Arm A (M)")
    p_avid.add_argument("--kd-b", type=float, default=1e-8, help="KD for Arm B (M)")
    p_avid.add_argument("--linker-len", type=int, default=15, help="Linker length (amino acids)")
    p_avid.add_argument("--c-eff", type=float, default=None, help="Explicit C_eff override (M)")
    p_avid.add_argument("--json", action="store_true", help="Output JSON")

    # Subcommand: tcell-lysis
    p_lysis = subparsers.add_parser("tcell-lysis", help="Simulate BiTE T-cell mediated cell lysis")
    p_lysis.add_argument("--construct-id", type=str, default="BiTE-001", help="Construct identifier")
    p_lysis.add_argument("--name", type=str, default="Blinatumomab-like", help="Construct name")
    p_lysis.add_argument("--tumor-antigen", type=str, default="CD19", help="Tumor target antigen")
    p_lysis.add_argument("--tumor-kd", type=float, default=1e-9, help="Tumor antigen KD (M)")
    p_lysis.add_argument("--cd3-kd", type=float, default=1e-7, help="CD3 effector arm KD (M)")
    p_lysis.add_argument("--tumor-density", type=float, default=50000.0, help="Receptors per tumor cell")
    p_lysis.add_argument("--cd3-density", type=float, default=40000.0, help="CD3 copies per T cell")
    p_lysis.add_argument("--et-ratio", type=float, default=5.0, help="Effector-to-target ratio")
    p_lysis.add_argument("--ab-conc", type=float, default=1e-10, help="Antibody concentration (M)")
    p_lysis.add_argument("--hours", type=float, default=24.0, help="Incubation time (hours)")
    p_lysis.add_argument("--alpha", type=float, default=1.0, help="Cooperativity factor")
    p_lysis.add_argument("--json", action="store_true", help="Output JSON")

    # Subcommand: fit-spr
    p_spr = subparsers.add_parser("fit-spr", help="Run SPR global 1:1 multi-cycle curve fitting")
    p_spr.add_argument("--json", action="store_true", help="Output JSON")

    # Subcommand: optimize-linker
    p_link = subparsers.add_parser("optimize-linker", help="Optimize linker length for epitope geometry")
    p_link.add_argument("--distance", type=float, default=4.5, help="Target epitope distance (nm)")
    p_link.add_argument("--kd-a", type=float, default=1e-9, help="KD for Arm A (M)")
    p_link.add_argument("--kd-b", type=float, default=1e-8, help="KD for Arm B (M)")
    p_link.add_argument("--json", action="store_true", help="Output JSON")

    # Subcommand: interactive
    subparsers.add_parser("interactive", help="Launch interactive analytics REPL")

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "ternary":
        return cmd_ternary(args)
    elif args.command == "titration":
        return cmd_titration(args)
    elif args.command == "avidity":
        return cmd_avidity(args)
    elif args.command == "tcell-lysis":
        return cmd_tcell_lysis(args)
    elif args.command == "fit-spr":
        return cmd_fit_spr(args)
    elif args.command == "optimize-linker":
        return cmd_optimize_linker(args)
    elif args.command == "interactive":
        return cmd_interactive()
    return 0


if __name__ == "__main__":
    sys.exit(main())
