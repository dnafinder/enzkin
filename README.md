⚗️ enzkin — Enzyme Kinetics Analysis in MATLAB

enzkin.m is a MATLAB function for estimating Michaelis–Menten kinetic parameters (Km and Vmax) from initial-rate enzyme kinetics data.
It combines classical linearized plots (Lineweaver–Burk, Hanes–Woolf, Eadie–Hofstee, Scatchard) with non-linear regression fits to provide a comprehensive view of the enzyme’s kinetic behaviour.

✨ Features

Accepts experimental pairs of substrate concentration [S] and initial velocity v

Performs 6 different regression approaches:

Lineweaver–Burk (1/v vs 1/[S])

Hanes–Woolf ([S]/v vs [S])

Eadie–Hofstee (v vs v/[S])

Scatchard (v/[S] vs v)

Logarithmic non-linear Michaelis–Menten fit

Hyperbolic (standard) non-linear Michaelis–Menten fit

Computes Km and Vmax with:

point estimate

standard error

lower and upper confidence bounds

Prints detailed regression summaries to the MATLAB Command Window

Produces a multi-panel figure with:

Michaelis–Menten hyperbolic plot

Lineweaver–Burk plot

Hanes–Woolf plot

log(Michaelis–Menten) plot

Eadie–Hofstee plot

Scatchard plot

🧬 Background

In enzyme kinetics, the Michaelis–Menten model relates the initial rate v to substrate concentration [S] via:

v = (Vmax · [S]) / (Km + [S])

where:

Vmax = maximum rate (asymptote)

Km = Michaelis–Menten constant (substrate concentration at v = Vmax/2)

While accurate parameter estimation should ideally rely on non-linear regression, classical linear transformations (Lineweaver–Burk, Hanes–Woolf, etc.) are still widely used for teaching and for quick diagnostic checks.
enzkin.m brings both worlds together in a single workflow.

🚀 Usage (MATLAB)

Basic example:

Define substrate concentrations and initial rates as row vectors:
S = [S₁ S₂ … Sₙ]
v = [v₁ v₂ … vₙ]

Call the function:
enzkinout = enzkin(S, v)

Output structure:

enzkinout.KM → 6 × 4 matrix with Km estimates:
[value, standard error, lower CI, upper CI]

enzkinout.VMAX → 6 × 4 matrix with Vmax estimates:
in the same order (linearisations + non-linear fits)

The function also prints:

slope, intercept, R, p-values, confidence intervals

summary tables for all 6 methods

And opens a figure window with 6 subplots (Michaelis–Menten and all linearised plots).

📦 Requirements

MATLAB

Curve Fitting Toolbox (for fit, fittype)

The custom regression function myregr by the same author, available at:
https://github.com/dnafinder/myregr

Make sure myregr.m is on your MATLAB path before calling enzkin.

📚 Citation

If you use this function in teaching, research, or publications, please cite:

Cardillo G. (2010). enzkin.m – A tool to estimate Michaelis–Menten kinetic parameters using multiple linear and non-linear regressions.
GitHub: https://github.com/dnafinder/enzkin

🔑 License

See the LICENSE file in this repository for licensing details.

👤 Author

Giuseppe Cardillo
Email: giuseppe.cardillo.75@gmail.com
