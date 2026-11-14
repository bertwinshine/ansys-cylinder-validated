#!/usr/bin/env python3
"""
Merge the ANSYS exports of the thick-walled pressure vessel, compare hoop and
radial stress with the Lame solution, and write pressure_vessel_reference.csv.

Inner and outer radius are read from the node coordinates.

Usage:
    python merge_and_validate.py
    python merge_and_validate.py --exports ./ansys_exports --pressure 3e5

Expected files in the exports directory (ANSYS 'Node Number<TAB>Value'):
    LOCX.txt, LOCY.txt, LOCZ.txt          undeformed coordinates (m)
    hoop_stress.txt                       Normal stress, cylindrical Y = hoop  (Pa)
    Radial_stress.txt                     Normal stress, cylindrical X = radial (Pa)
    Equivalent_stress.txt                 von Mises (Pa)
    Total_deformation.txt                  total deformation (m)

Bertwin Kurisinkal Shine, TU Chemnitz. Written with help from an AI assistant.
"""

import argparse
import os
import numpy as np
import pandas as pd


# Lame solution for a thick-walled cylinder
def lame_stresses(r, a, b, p_i, p_o=0.0):
    """Exact radial and hoop stress for a thick-walled cylinder.

    r : radius/radii (a <= r <= b);  a, b : inner/outer radius
    p_i, p_o : internal / external pressure
    """
    r = np.asarray(r, float)
    a2, b2, r2 = a * a, b * b, r * r
    denom = b2 - a2
    sr = (p_i * a2 - p_o * b2) / denom - (p_i - p_o) * a2 * b2 / (denom * r2)
    st = (p_i * a2 - p_o * b2) / denom + (p_i - p_o) * a2 * b2 / (denom * r2)
    return sr, st


def load(exports_dir, name, valcol):
    """Read one ANSYS export -> DataFrame[node, valcol]."""
    path = os.path.join(exports_dir, name)
    df = pd.read_csv(path, sep="\t", engine="python")
    df.columns = ["node", valcol]
    df[valcol] = pd.to_numeric(df[valcol], errors="coerce")
    return df


def main():
    ap = argparse.ArgumentParser(description="Validate ANSYS pressure-vessel "
                                             "result against the Lame solution.")
    ap.add_argument("--exports", default="ansys_exports",
                    help="folder containing the ANSYS .txt exports "
                         "(default: ansys_exports)")
    ap.add_argument("--pressure", type=float, default=3.0e5,
                    help="internal pressure in Pa (default: 3e5)")
    ap.add_argument("--out", default="pressure_vessel_reference.csv",
                    help="output CSV path")
    ap.add_argument("--tol", type=float, default=2.0,
                    help="acceptance tolerance, %% of peak hoop stress (default 2)")
    args = ap.parse_args()

    P_IN = args.pressure
    ex = args.exports

    # ---- load + merge ----------------------------------------------------
    x = load(ex, "LOCX.txt", "x")
    y = load(ex, "LOCY.txt", "y")
    z = load(ex, "LOCZ.txt", "z")
    hp = load(ex, "hoop_stress.txt", "hoop")
    rd = load(ex, "Radial_stress.txt", "radial")
    vm = load(ex, "Equivalent_stress.txt", "von_mises")
    df_def = load(ex, "Total_deformation.txt", "deformation")

    df = (x.merge(y, on="node").merge(z, on="node")
            .merge(rd, on="node").merge(hp, on="node")
            .merge(vm, on="node").merge(df_def, on="node"))
    df["radius"] = np.sqrt(df["x"] ** 2 + df["y"] ** 2)

    # ---- auto-detect geometry -------------------------------------------
    a_in = round(df["radius"].min(), 4)
    b_out = round(df["radius"].max(), 4)
    print(f"Merged {len(df):,} nodes.")
    print(f"Auto-detected geometry: inner a = {a_in} m, outer b = {b_out} m, "
          f"thickness = {b_out - a_in:.4f} m")
    print(f"Axial length (z): {df['z'].min():.3f} .. {df['z'].max():.3f} m")
    print(f"Internal pressure: {P_IN:.3e} Pa\n")

    # ---- validate (central 80% of length, skip end-face edge effects) ---
    zmin, zmax = df["z"].min(), df["z"].max()
    zlo = zmin + 0.10 * (zmax - zmin)
    zhi = zmax - 0.10 * (zmax - zmin)
    core = df[(df["z"] >= zlo) & (df["z"] <= zhi)].copy()

    sr_ref, st_ref = lame_stresses(core["radius"].values, a_in, b_out, P_IN)
    scale = np.max(np.abs(st_ref))
    err_hoop = np.abs(core["hoop"].values - st_ref) / scale * 100
    err_radial = np.abs(core["radial"].values - sr_ref) / scale * 100

    print("=" * 64)
    print("VALIDATION AGAINST EXACT LAME SOLUTION  (central 80% of length)")
    print("=" * 64)
    print(f"  nodes checked     : {len(core):,}")
    print(f"  mean hoop  error  : {err_hoop.mean():.3f} %")
    print(f"  max  hoop  error  : {err_hoop.max():.3f} %")
    print(f"  mean radial error : {err_radial.mean():.3f} %")
    print(f"  max  radial error : {err_radial.max():.3f} %")
    passed = max(err_hoop.mean(), err_radial.mean()) <= args.tol
    print(f"  tolerance (mean)  : {args.tol:.1f} %")
    print(f"  RESULT            : {'PASS' if passed else 'FAIL'}")
    print()

    # ---- wall spot-checks ------------------------------------------------
    inner = df[np.isclose(df["radius"], a_in, atol=3e-3)]
    outer = df[np.isclose(df["radius"], b_out, atol=3e-3)]
    sr_in, st_in = lame_stresses(a_in, a_in, b_out, P_IN)
    sr_out, st_out = lame_stresses(b_out, a_in, b_out, P_IN)
    print("  wall-value spot check (ANSYS mean  vs  Lame exact):")
    print(f"    hoop  inner : {inner['hoop'].mean():.4e}  vs {st_in:.4e} Pa")
    print(f"    hoop  outer : {outer['hoop'].mean():.4e}  vs {st_out:.4e} Pa")
    print(f"    radial inner: {inner['radial'].mean():.4e}  vs {sr_in:.4e} Pa (= -p)")
    print(f"    radial outer: {outer['radial'].mean():.4e}  vs {sr_out:.4e} Pa")
    print(f"    von Mises max (field): {df['von_mises'].max():.4e} Pa")
    print()

    # ---- write deliverable CSV ------------------------------------------
    out_cols = ["node", "x", "y", "z", "radius",
                "radial", "hoop", "von_mises", "deformation"]
    df_out = df[out_cols].sort_values("node").reset_index(drop=True)
    df_out.to_csv(args.out, index=False)
    print(f"Wrote certified reference field -> {args.out}")
    print(f"  rows: {len(df_out):,}   columns: {', '.join(out_cols)}")


if __name__ == "__main__":
    main()
