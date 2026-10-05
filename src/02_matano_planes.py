"""Determine Matano planes and draw one combined shaded figure."""
from workflow import paths, load_profiles, matano_planes

if __name__ == "__main__":
    excel, output = paths(__doc__)
    matano_planes(load_profiles(excel), output)
