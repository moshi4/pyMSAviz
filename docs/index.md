# pyMSAviz

![Python3](https://img.shields.io/badge/Language-Python3-steelblue)
![OS](https://img.shields.io/badge/OS-_Windows_|_Mac_|_Linux-steelblue)
![License](https://img.shields.io/badge/License-MIT-steelblue)
[![Latest PyPI version](https://img.shields.io/pypi/v/pymsaviz.svg)](https://pypi.python.org/pypi/pymsaviz)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/pymsaviz.svg?color=green)](https://anaconda.org/bioconda/pymsaviz)  

## Overview

pyMSAviz is a MSA(Multiple Sequence Alignment) visualization python package for sequence analysis implemented based on matplotlib.
This package is developed for the purpose of easily and beautifully plotting MSA in Python.
It also implements the functionality to add markers, text annotations, highlights to specific positions and ranges in MSA.
pyMSAviz was developed inspired by [Jalview](https://www.jalview.org/) and [ggmsa](https://github.com/YuLab-SMU/ggmsa).

<figure markdown>
  ![example.png](./images/api_example01.png)
  <figcaption>Fig.1 Simple visualization result</figcaption>
</figure>

<figure markdown>
  ![example.png](./images/api_example03.png)
  <figcaption>Fig.2 Customized visualization result</figcaption>
</figure>

## Installation

`Python 3.9 or later` is required for installation.

**Install PyPI package:**

    pip install pymsaviz

**Install bioconda package:**

    conda install -c conda-forge -c bioconda pymsaviz
