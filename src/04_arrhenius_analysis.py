"""Fit mean coefficients versus reciprocal absolute temperature."""
from workflow import paths, load_profiles, arrhenius

if __name__ == "__main__":
    excel, output = paths(__doc__)
    arrhenius(load_profiles(excel), output)
