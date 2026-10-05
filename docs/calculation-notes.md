# Calculation notes

The main README follows the concentration-profile, Matano-plane, coefficient, and Arrhenius sequence in the supplied diffusion report. Only its relevant theory (sections 2.3–2.5) and Co–Al interdiffusion discussion (section 4.2 and its Arrhenius subsection) were used. The original Matano reference is listed in the README. Unrelated Cu–Zn experiments are outside this repository's scope.

## Current calculation versus report

The source report gives Q = 317.93 kJ/mol and D₀ = 5.988 × 10⁻³ m²/s. Regression of its supplied mean coefficients reproduces those parameters. However, the current unsmoothed position-space calculation from the uploaded workbook gives Q = 334.75 kJ/mol and D₀ = 2.56963 × 10⁻² m²/s. The report's numbers have not been substituted into recalculated figures or tables.

| Temperature | Earlier reported mean (m²/s) | Current raw-profile mean (m²/s) |
|---|---:|---:|
| 1100 °C | 4.68966708e−15 | 4.51026338e−15 |
| 1200 °C | 3.37629221e−14 | 3.85468597e−14 |
| 1300 °C | 1.61098946e−13 | 1.86222651e−13 |

The supplied example image guides the combined Matano figure's layout. Its labels reflect the earlier calculation (501.10, 1010.31, 789.38 μm). The regenerated plot uses the current mass-balance values (501.11, 1010.43, 789.40 μm). Dashed-line colors identify the corresponding profile, and overlapping shading is illustrative rather than a numerical area test.

## Integration method

The original notebook attempted inverse interpolation x(C) on noisy concentration data. These profiles contain repeated concentrations and local reversals; reversing arrays does not make the inverse single-valued. The current default integrates in position space using integration by parts. This preserves the BM mass-balance relation and the raw observations without smoothing or sorting by concentration.

The earlier inverse function remains in `src/diffusion.py` for numerical comparison on suitable monotonic inputs, with a guard that rejects nonmonotonic inputs. Stage scripts use only the position-space function.

Positions are metres internally, times are seconds, and coefficients are m²/s. The raw derivative uses NumPy finite differences. The inherited very small gradient threshold (10⁻⁶ wt%/m) is a division safeguard, not an experimental noise threshold. The mean is an arithmetic mean of samples in 0.5–2.5 wt% Al. All in-range estimates are finite and positive in this run; complete exported arrays retain other invalid or negative values.

## Workbook formulas

The workbook is retained unchanged. Its raw columns A and B and annealing-time parameters are read; existing coefficient formulas are not used.

- Each temperature sheet's L5 uses M5/M6, a concentration-weighted centroid expression instead of the Matano mass balance. Using its cumulative integral E, the equivalent mass-balance formula would be `=A1001-E1001/(B1001-B2)` for 1100C and 1200C, and `=A801-E801/(B801-B2)` for 1300C.
- `Arrhenius!F4` uses the intercept and activation-energy cells instead of slope and intercept. It should be `=$B$8*C4+$B$9`, filled through F6. This is an identified source issue, not an edit made to the supplied workbook.

## Interpretation

Measured endpoints are used as terminal concentrations. Profile-end fluctuations and incomplete plateaus can affect the calculation. The use of wt% directly assumes an appropriate constant-density concentration balance. No endpoint-sensitivity study, density correction, or experimental uncertainty model is included.

The calculated coefficients do not increase monotonically with composition at every temperature. Coordinate offsets do not alone demonstrate a moving physical interface, and three approximately collinear Arrhenius points do not identify a unique diffusion mechanism. These stronger claims from the report were not carried into the README.

## Verification

Numerical tests use an analytical constant-diffusivity solution, translated coordinates, reversed concentration orientation, known Arrhenius parameters, and invalid inputs. The full workflow and individual stage entry points are checked against the same supplied workbook. Visual inspection covers the four output plots. These checks establish computational consistency under the stated assumptions, not independent experimental validation.
