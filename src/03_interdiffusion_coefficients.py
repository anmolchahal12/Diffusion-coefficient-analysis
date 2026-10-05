"""Calculate coefficients and the mean table."""
from workflow import paths, load_profiles, interdiffusion

if __name__ == "__main__":
    excel, output = paths(__doc__)
    interdiffusion(load_profiles(excel), output)
