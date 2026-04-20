# Resonant Vibrations of Viral Capsids: ENM Code

(https://github.com/EdilsonSilva23/viral-capsid-resonances)

This repository contains the Python implementation of the Elastic Network Model (ENM) described in:

**Silva, E.** *Resonant Vibrations of Viral Capsids: A Unified Framework from Lamb Theory to Genome-Induced Prestress*. Submitted to *Soft Matter* (Royal Society of Chemistry).

## Overview

The code constructs an ENM from a viral capsid PDB structure, calibrates the spring constant using continuum elasticity (exact thin‑shell denominator), and computes the vibrational spectrum. It reproduces the symmetry decomposition and frequency predictions reported in the manuscript.

## Requirements

- Python 3.10+
- NumPy
- SciPy
- BioPython
- Matplotlib

Install dependencies with:

```bash
pip install numpy scipy biopython matplotlib