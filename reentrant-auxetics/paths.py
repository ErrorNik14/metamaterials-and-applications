"""Folder locations shared by every script in this folder."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
INPUTS = ROOT / "inputs"                                   # STEP files + Fukuchi data
OUTPUTS = ROOT / "outputs"                                 # everything the scripts write
ANSYS_RESULTS = ROOT / "ansys" / "results"    # hand-run ANSYS CSVs

OUTPUTS.mkdir(exist_ok=True)
