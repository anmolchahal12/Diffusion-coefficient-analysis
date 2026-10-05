"""Plot concentration versus position at all temperatures."""
from workflow import paths, load_profiles, concentration_profiles

if __name__ == "__main__":
    excel, output = paths(__doc__)
    concentration_profiles(load_profiles(excel), output)
