"""Workbook loading, calculation, and export helpers shared by the stage scripts."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from diffusion import prepare_profile, position_space_diffusivity, mean_in_interval, fit_arrhenius

ROOT = Path(__file__).resolve().parents[1]
CONDITIONS = [(1100, 201.0), (1200, 151.0), (1300, 76.0)]
COLORS = ["#1976b9", "#e68418", "#299441"]

def paths(description):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--excel", type=Path, default=ROOT / "data/Boltzmann_Matano_CoAl.xlsx")
    parser.add_argument("--output", type=Path, default=ROOT / "results")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    return args.excel, args.output

def load_profiles(excel):
    """Read positions in μm and concentrations in wt% Al; validate annealing times."""
    profiles = {}
    with pd.ExcelFile(excel) as book:
        for temperature, hours in CONDITIONS:
            frame = pd.read_excel(book, sheet_name=f"{temperature}C")
            if list(frame.columns[:2]) != ["x (μm)", "Wt% Al"]:
                raise ValueError("Expected position in μm and concentration in wt% Al.")
            if frame.iloc[1, 10] != "Anneal time (h)" or float(frame.iloc[1, 11]) != hours:
                raise ValueError(f"Annealing time differs from configured {hours} h.")
            x = pd.to_numeric(frame.iloc[:, 0], errors="raise").to_numpy(float)*1e-6
            c = pd.to_numeric(frame.iloc[:, 1], errors="raise").to_numpy(float)
            x, c = prepare_profile(x, c)
            profiles[temperature] = {"x_m": x, "C": c, "hours": hours}
    return profiles

def save_figure(fig, output, name):
    fig.savefig(output / f"{name}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

def concentration_profiles(profiles, output):
    """Stage 1: raw concentration versus measured position."""
    fig, ax = plt.subplots(figsize=(9, 5.8), constrained_layout=True)
    for (temperature, p), color in zip(profiles.items(), COLORS):
        ax.plot(p["x_m"]*1e6, p["C"], color=color, label=f"{temperature} °C", lw=1.7)
    ax.set(xlabel="Position (μm)", ylabel="Aluminum concentration (wt% Al)", title="Co–Al concentration profiles")
    ax.legend()
    ax.grid(alpha=0.22)
    save_figure(fig, output, "01_concentration_profiles")

def matano_planes(profiles, output):
    """Stage 2: one combined figure; raw-profile shading is illustrative."""
    from matplotlib.patches import Patch
    fig, ax = plt.subplots(figsize=(11, 7), constrained_layout=True)
    rows, handles = [], []
    for (temperature, p), color in zip(profiles.items(), COLORS):
        x, c = p["x_m"], p["C"]
        xm = float(x[-1] - np.trapezoid(c-c[0], x=x)/(c[-1]-c[0]))
        xu, xmu = x*1e6, xm*1e6
        line, = ax.plot(xu, c, color=color, lw=1.8, label=f"{temperature} °C profile", zorder=3)
        plane = ax.axvline(xmu, color=color, ls="--", lw=1.4, label=f"Matano {temperature} °C: {xmu:.2f} μm", zorder=4)
        ax.fill_betweenx(c, xu, xmu, where=xu<=xmu, color="#88cce5", alpha=0.16, interpolate=True)
        ax.fill_betweenx(c, xu, xmu, where=xu>=xmu, color="#ee9999", alpha=0.16, interpolate=True)
        handles.extend([line, plane])
        rows.append({"temperature_C": temperature, "annealing_hours": p["hours"], "matano_position_um": xmu})
    handles.extend([Patch(facecolor="#88cce5", alpha=0.4, label="Left of each Matano plane"),
                    Patch(facecolor="#ee9999", alpha=0.4, label="Right of each Matano plane")])
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=9, frameon=False)
    ax.set(xlabel="Position (μm)", ylabel="Aluminum concentration (wt% Al)", title="Concentration profiles with Matano planes")
    ax.grid(alpha=0.22)
    save_figure(fig, output, "02_matano_planes")
    pd.DataFrame(rows).to_csv(output / "matano_positions.csv", index=False)

def calculate_coefficients(profiles):
    coefficients, means, diagnostics = {}, [], {}
    for temperature, p in profiles.items():
        x, c, d, xm, diag = position_space_diffusivity(p["x_m"], p["C"], p["hours"]*3600)
        mean, count = mean_in_interval(c, d)
        coefficients[temperature] = (x, c, d)
        diagnostics[str(temperature)] = diag
        means.append({"temperature_C": temperature, "annealing_hours": p["hours"],
                      "mean_D_m2_s": mean, "n_samples": count, "matano_position_um": xm*1e6})
    return coefficients, pd.DataFrame(means), diagnostics

def interdiffusion(profiles, output):
    """Stage 3: full D tables, interval plot, and sample arithmetic means."""
    coefficients, means, diagnostics = calculate_coefficients(profiles)
    fig, ax = plt.subplots(figsize=(9, 5.8), constrained_layout=True)
    for (temperature, (x, c, d)), color in zip(coefficients.items(), COLORS):
        pd.DataFrame({"position_m": x, "Al_wt_percent": c, "D_m2_s": d}).to_csv(output / f"diffusivity_{temperature}C.csv", index=False)
        selected = (c >= 0.5) & (c <= 2.5)
        shown = np.where(np.isfinite(d) & (d > 0), d, np.nan)
        ax.plot(c[selected], shown[selected], color=color, lw=1.5, label=f"{temperature} °C")
    ax.set(xlabel="Aluminum concentration (wt% Al)", ylabel="Interdiffusion coefficient (m²/s)",
           title="Composition-dependent interdiffusion coefficients", yscale="log")
    ax.legend()
    ax.grid(alpha=0.22, which="both")
    save_figure(fig, output, "03_interdiffusion_coefficients")
    means.to_csv(output / "mean_coefficients.csv", index=False)
    (output / "diagnostics.json").write_text(json.dumps(diagnostics, indent=2)+"\n")
    return means

def arrhenius(profiles, output, means=None):
    """Stage 4: recompute means when run independently, then fit against 1/T."""
    if means is None:
        _, means, _ = calculate_coefficients(profiles)
    means.to_csv(output / "mean_coefficients.csv", index=False)
    temp, d = means.temperature_C.to_numpy(), means.mean_D_m2_s.to_numpy()
    fit = fit_arrhenius(temp, d)
    inverse_t = 1/(temp+273.15)
    line = np.linspace(inverse_t.min(), inverse_t.max(), 100)
    fig, ax = plt.subplots(figsize=(9, 5.8), constrained_layout=True)
    ax.scatter(1000*inverse_t, np.log(d), color=COLORS, s=55, zorder=3, label="Mean coefficients")
    ax.plot(1000*line, fit["intercept"]+fit["slope_K"]*line, "--", color="#444444", label="Linear fit")
    for tc, tx, dy in zip(temp, inverse_t, np.log(d)):
        ax.annotate(f"{tc} °C", (1000*tx,dy), xytext=(8,8), textcoords="offset points", fontsize=9)
    ax.text(0.04, 0.08, f"Q = {fit['Q_kJ_mol']:.2f} kJ/mol\nD₀ = {fit['D0_m2_s']:.3e} m²/s\nR² = {fit['r_squared']:.6f}",
            transform=ax.transAxes, fontsize=10, bbox={"facecolor":"white", "edgecolor":"#dddddd", "pad":8})
    ax.set(xlabel="1000/T (K⁻¹)", ylabel="ln[D̄ / (1 m² s⁻¹)]", title="Arrhenius analysis of Co–Al interdiffusion")
    ax.margins(x=0.12, y=0.12)
    ax.legend(loc="upper right")
    ax.grid(alpha=0.22)
    save_figure(fig, output, "04_arrhenius")
    (output / "arrhenius_fit.json").write_text(json.dumps(fit, indent=2)+"\n")
    print(f"Q = {fit['Q_kJ_mol']:.2f} kJ/mol; D0 = {fit['D0_m2_s']:.6g} m²/s")
