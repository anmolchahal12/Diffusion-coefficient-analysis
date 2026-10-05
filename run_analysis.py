"""Run all four Co–Al analysis stages from the included workbook."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from workflow import paths, load_profiles, concentration_profiles, matano_planes, interdiffusion, arrhenius

if __name__ == "__main__":
    excel, output = paths(__doc__)
    profiles = load_profiles(excel)
    concentration_profiles(profiles, output)
    matano_planes(profiles, output)
    means = interdiffusion(profiles, output)
    arrhenius(profiles, output, means)
    print(f"Saved figures and tables to {output}")
