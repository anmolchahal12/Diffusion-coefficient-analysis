"""Shared Boltzmann–Matano calculations. Inputs: metres, wt% Al, seconds."""
import numpy as np
from scipy.integrate import simpson, cumulative_trapezoid
from scipy.interpolate import interp1d
from scipy.stats import linregress
R = 8.314
LOW, HIGH = 0.5, 2.5

def prepare_profile(position_m, concentration):
    """Return finite data ordered by position; never discard rows silently."""
    x, c = np.asarray(position_m, dtype=float), np.asarray(concentration, dtype=float)
    if x.ndim != 1 or c.ndim != 1 or x.shape != c.shape or len(x) < 5:
        raise ValueError("Need matching one-dimensional arrays with at least five points.")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(c)):
        raise ValueError("Non-numeric, missing, or infinite input: inspect the worksheet rows.")
    order = np.argsort(x)
    x, c = x[order], c[order]
    if np.any(np.diff(x) <= 0):
        raise ValueError("Duplicate position values: resolve them before differentiation.")
    return x, c


def boltzmann_matano_diffusivity(position_m, concentration, time_seconds):
    """Return all oriented samples, D (including NaNs), x_M, and diagnostics.

    Inputs: metres, wt% Al, seconds. Output D: m²/s. The constant-density/
    volume and concentration-basis assumptions still require experimental review.
    Simpson integration and raw finite-difference gradients follow the notebook.
    """
    if not np.isfinite(time_seconds) or time_seconds <= 0:
        raise ValueError("Annealing time must be finite and positive, in seconds.")
    x, c = prepare_profile(position_m, concentration)
    if c[0] > c[-1]:
        x, c = x[::-1], c[::-1]
    if np.any(np.diff(c) <= 0):
        raise ValueError(
            "Raw concentrations are not strictly monotonic (or contain repeats). "
            "The notebook's inverse x(C) is not well-defined for this profile. "
            "No smoothing or concentration sorting was applied. Review the raw "
            "profile and choose a documented treatment before calculating D."
        )
    x_m = float(simpson(x, x=c) / (c[-1] - c[0]))
    shifted = x - x_m
    gradient = np.gradient(c, shifted)
    # Inherited numerical threshold, in wt% Al per metre; not a noise criterion.
    valid_gradient = np.isfinite(gradient) & (np.abs(gradient) > 1e-6)
    if valid_gradient.sum() < 2:
        raise ValueError("Fewer than two usable gradients.")
    # No extrapolation: insufficient integral support is explicitly invalidated.
    inverse = interp1d(c[valid_gradient], shifted[valid_gradient],
                      kind="linear", bounds_error=False, fill_value=np.nan)
    integrals = np.full(c.shape, np.nan)
    for i in range(len(c)):
        if not valid_gradient[i]:
            continue
        if i == 0:
            integrals[i] = 0.0
        else:
            integrals[i] = simpson(inverse(c[:i + 1]), x=c[:i + 1])
    d = np.full(c.shape, np.nan)
    np.divide(-integrals, 2 * time_seconds * gradient,
              out=d, where=valid_gradient)
    selected = (c >= LOW) & (c <= HIGH)
    diagnostics = {
        "points_total": int(len(c)),
        "points_in_interval": int(selected.sum()),
        "nonfinite_D_total": int((~np.isfinite(d)).sum()),
        "nonpositive_D_total": int((np.isfinite(d) & (d <= 0)).sum()),
        "invalid_D_in_interval": int((selected & (~np.isfinite(d) | (d <= 0))).sum()),
        "matano_position_um": x_m * 1e6,
        "smoothing": "none",
    }
    return x, c, d, x_m, diagnostics


def mean_in_interval(c, d):
    """Sample arithmetic mean, matching the original notebook's convention."""
    if c.min() > LOW or c.max() < HIGH:
        raise ValueError("Profile does not span the entire 0.5–2.5 wt% Al interval.")
    selected = (c >= LOW) & (c <= HIGH)
    values = d[selected]
    if len(values) < 2 or not np.all(np.isfinite(values) & (values > 0)):
        raise ValueError("Invalid D values in the averaging interval; no positivity-filtered mean produced.")
    return float(np.mean(values)), int(len(values))


def position_space_diffusivity(position_m, concentration, time_seconds):
    """Raw BM estimate without constructing the potentially multivalued x(C).

    Let a = C-C_left and A(x) = integral(a dx). Integration by parts gives
    x_M = x_right - A(x_right)/(C_right-C_left) and
    integral((x-x_M) dC) = (x-x_M)*a - A(x).
    Trapezoid integration follows position order, preserving concentration noise.
    No curve fitting, smoothing, concentration sorting, or clipping is performed.
    Results assume representative terminal compositions and the stated
    constant-density/concentration formulation. They remain raw estimates.
    """
    if not np.isfinite(time_seconds) or time_seconds <= 0:
        raise ValueError("Annealing time must be finite and positive, in seconds.")
    x, c = prepare_profile(position_m, concentration)
    if c[-1] == c[0]:
        raise ValueError("Terminal concentrations are equal; the Matano plane is undefined.")
    a = c - c[0]
    area = cumulative_trapezoid(a, x=x, initial=0)
    xm = float(x[-1] - area[-1] / a[-1])
    integral = (x-xm)*a - area
    gradient = np.gradient(c, x)
    valid = np.isfinite(gradient) & (np.abs(gradient) > 1e-6)
    d = np.full(c.shape, np.nan)
    np.divide(-integral, 2*time_seconds*gradient, out=d, where=valid)
    selected = (c >= LOW) & (c <= HIGH)
    direction = np.sign(c[-1]-c[0])
    diagnostics = {
        "points_total": int(len(c)),
        "points_in_interval": int(selected.sum()),
        "nonfinite_D_total": int((~np.isfinite(d)).sum()),
        "nonpositive_D_total": int((np.isfinite(d) & (d <= 0)).sum()),
        "invalid_D_in_interval": int((selected & (~np.isfinite(d) | (d <= 0))).sum()),
        "local_steps_against_overall_trend": int((np.diff(c)*direction < 0).sum()),
        "repeated_concentration_values": int(len(c)-len(np.unique(c))),
        "matano_position_um": xm*1e6,
        "end_integral_m_wt_percent": float(integral[-1]),
        "concentration_at_left_wt_percent": float(c[0]),
        "concentration_at_right_wt_percent": float(c[-1]),
        "smoothing": "none",
        "method": "position-space integration by parts; trapezoid quadrature",
    }
    return x, c, d, xm, diagnostics


def fit_arrhenius(temperatures_c, mean_d):
    """Fit ln(D / [1 m²/s]) versus 1/T, preserving R=8.314 from the notebook."""
    temp, d = np.asarray(temperatures_c, float), np.asarray(mean_d, float)
    if (temp.ndim != 1 or temp.shape != d.shape or len(temp) < 3
            or len(np.unique(temp)) != len(temp)
            or not np.all(np.isfinite(temp) & (temp > -273.15))
            or not np.all(np.isfinite(d) & (d > 0))):
        raise ValueError("Need at least three distinct finite temperatures and positive finite means.")
    fit = linregress(1 / (temp + 273.15), np.log(d))
    return {"Q_kJ_mol": float(-fit.slope * R / 1000),
            "D0_m2_s": float(np.exp(fit.intercept)),
            "slope_K": float(fit.slope), "intercept": float(fit.intercept),
            "r_squared": float(fit.rvalue ** 2), "gas_constant_J_mol_K": R}


