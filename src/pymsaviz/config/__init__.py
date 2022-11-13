from __future__ import annotations

import csv
from enum import IntEnum, auto
from pathlib import Path

###########################################################
# Color Schemes Config
###########################################################


def get_color_schemes() -> dict[str, dict[str, str]]:
    """Get color schemes

    Returns
    -------
    name2color_scheme : dict[str, dict[str, str]]
        Color schemes dict
    """
    COLOR_SCHEMES_FILE = Path(__file__).parent / "color_schemes.tsv"
    name2color_scheme = {}
    with open(COLOR_SCHEMES_FILE) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        letters = header[1:]
        for row in reader:
            name, colors = row[0], row[1:]
            color_scheme = {}
            for letter, color in zip(letters, colors):
                color_scheme[letter] = color
            name2color_scheme[name] = color_scheme
    return name2color_scheme


COLOR_SCHEMES = get_color_schemes()

###########################################################
# Plot Config
###########################################################


class AxesType(IntEnum):
    """Plot axes type enum"""

    MSA = auto()
    CONSENSUS = auto()
    SPACE = auto()
    WRAP_SPACE = auto()


###########################################################
# Example MSA Dataset
###########################################################


def get_msa_testdata(name: str = "MRGPRG.fa") -> Path:
    """Get MSA testdata file

    List of MSA testdata filename
    - `HIGD2A.fa` (6 species genes, 118 alignment length)
    - `MRGPRG.fa` (6 species genes, 289 alignment length)

    Parameters
    ----------
    name : str, optional
        Testdata name

    Returns
    -------
    msa_testdata_file : Path
        MSA testdata file
    """
    testdata_dir = Path(__file__).parent / "testdata"
    dataset_files = testdata_dir.glob("*")
    name2dataset_file = {f.name: f for f in dataset_files}
    if name not in name2dataset_file:
        err_msg = f"Dataset name = '{name}' not found. "
        err_msg += f"Available testdata name = {list(name2dataset_file.keys())}"
        raise ValueError(err_msg)
    return name2dataset_file[name]
