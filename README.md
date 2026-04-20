# Resonant Vibrations of Viral Capsids: ENM Code

(https://github.com/EdilsonSilva23/viral-capsid-resonances)

This repository contains the Python implementation of the Elastic Network Model (ENM) described in:

**Silva, E.** *Resonant Vibrations of Viral Capsids: A Unified Framework from Lamb Theory to Genome-Induced Prestress*. *Soft Matter* (submitted).

## Overview

The code constructs an ENM from a viral capsid PDB structure, calibrates the spring constant using continuum elasticity (exact thin‑shell denominator), and computes the vibrational spectrum. It reproduces the symmetry decomposition and frequency predictions reported in the manuscript.

---

## How to Obtain and Run the Code

### 1. Clone or download the repository

```bash
git clone https://github.com/edilsonsilva/viral-capsid-resonances.git
cd viral-capsid-resonances
Alternatively, download the ZIP archive from the GitHub page or directly from the Zenodo record.

2. Set up the Python environment
The code requires Python 3.10+ and the following packages:

numpy

scipy

biopython

matplotlib

You can install the dependencies using pip and the provided requirements.txt file:

bash
pip install -r requirements.txt
3. Obtain the capsid structure (PDB file)
The example provided uses Cowpea Chlorotic Mottle Virus (CCMV), PDB ID 1CWP.
You can download the structure directly from the Protein Data Bank:

bash
wget https://files.rcsb.org/download/1CWP.pdb
Or visit https://www.rcsb.org/structure/1CWP and download the PDB Format file manually.

Note: The code works with any viral capsid PDB file containing Cα atoms. Replace 1CWP.pdb with the path to your file.

4. Run the ENM analysis
Execute the script with the PDB file as an argument:

bash
python enm_capsid.py 1CWP.pdb
Output:

Continuum frequency predictions (exact vs. approximate denominators).

Fractional prestress shifts for reference viruses.

Calibrated spring constant k_cal and effective node mass.

The lowest 10 non‑trivial ENM mode frequencies.

Supplementary Figure S1 (Figure_S1_ENM_CCMV.pdf / .png): vibrational spectrum showing icosahedral degeneracy patterns.

If no PDB file is supplied, the script only computes the continuum predictions (useful for quick parameter testing).

Key Equations Implemented
The code uses the exact thin‑shell inertial denominator:

text
ρ_denom = 3 ρ_shell h / R + 3 ρ_f
and the generalized resonance equation (main text Eq. 7):

text
f = (ξ / 2πR) √[ (E + γ/R) / (3ρ_shell h/R + 3ρ_f) ]
with Lamb eigenvalue ξ = 2.08 for the breathing mode (ν = 0.3).

Citation
If you use this code in your research, please cite the accompanying manuscript:

Silva, E. Resonant Vibrations of Viral Capsids: A Unified Framework from Lamb Theory to Genome-Induced Prestress. Soft Matter (submitted).

and the Zenodo DOI of this repository.

License
This code is released under the MIT License. See the LICENSE file for details.

Contact
Edilson Silva – edilsonsilva@hotmail.com
