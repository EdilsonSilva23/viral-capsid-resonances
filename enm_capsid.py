"""
Resonant Vibrations of Viral Capsids: A Unified Framework
Edilson Silva — Soft Matter submission

Elastic Network Model (ENM) for CCMV (PDB: 1CWP)

Exact version:
- Inertial denominator: 3*rho_shell*h/R + 3*rho_f (Supplementary Note 6)
- Lamb eigenvalue xi = 2.08 (nu = 0.3, breathing mode)
- Dynamical matrix D = K_scaled / m_node (units: rad^2/s^2)

Usage:
    python enm_capsid.py 1cwp.pdb
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigsh

# ============================================================
# PHYSICAL PARAMETERS
# ============================================================
xi_10 = 2.08          # Lamb eigenvalue, breathing mode (l=0), nu=0.3
xi_12 = 2.65          # Lamb eigenvalue, quadrupolar mode (l=2), nu=0.3
rho_shell = 1200.0    # kg/m^3, bulk protein shell density
rho_f = 1000.0        # kg/m^3, water density at 298 K
pi = np.pi

# ============================================================
# INERTIAL DENOMINATORS
# ============================================================
def exact_denominator(R_m, h_m, rho_s=rho_shell, rho_fluid=rho_f):
    """
    Exact thin-shell inertial denominator (kg/m^3).
    Eq. (S1): rho_denom = 3*rho_shell*h/R + 3*rho_f
    """
    return 3.0 * rho_s * h_m / R_m + 3.0 * rho_fluid

def approx_denominator(C_val=3.0, rho_s=rho_shell, rho_fluid=rho_f):
    """
    Parametrized form: rho_shell + C*rho_f.
    C ~ 3 for thin-shell breathing mode (geometry-dependent; Supplementary Note 6).
    """
    return rho_s + C_val * rho_fluid

# ============================================================
# FREQUENCY PREDICTIONS
# ============================================================
def predict_frequency_exact(R_m, E_Pa, h_m=2e-9, gamma_Nm=0.0,
                            rho_s=rho_shell, rho_fluid=rho_f, xi=xi_10):
    """Frequency using exact thin-shell denominator (main text Eq. 7)."""
    denom = exact_denominator(R_m, h_m, rho_s, rho_fluid)
    E_eff = E_Pa + gamma_Nm / R_m
    return (xi / (2.0 * pi * R_m)) * np.sqrt(E_eff / denom)

def predict_frequency_approx(R_m, E_Pa, C_val=3.0, gamma_Nm=0.0, xi=xi_10):
    """Frequency using parametrized form (conceptual comparison only)."""
    denom = approx_denominator(C_val)
    E_eff = E_Pa + gamma_Nm / R_m
    return (xi / (2.0 * pi * R_m)) * np.sqrt(E_eff / denom)

def prestress_shift(P_Pa, E_Pa):
    """Fractional frequency shift from genome pressure (first-order, Eq. 8)."""
    return P_Pa / (4.0 * E_Pa)

def propagate_error(sigma_R_frac, sigma_E_frac, sigma_rho_frac=0.05):
    """Fractional frequency uncertainty from error propagation (Eq. 9)."""
    return np.sqrt(sigma_R_frac**2 +
                   (0.5 * sigma_E_frac)**2 +
                   (0.5 * sigma_rho_frac)**2)

# ============================================================
# ENM CONSTRUCTION FROM PDB
# ============================================================
def build_enm(pdb_file, cutoff_angstrom=10.0):
    """
    Parse PDB, extract CA atoms, build elastic network model.
    Coordinates converted from Angstrom to metres.
    Returns (coords_m, K_csr, N) where K is the dimensionless
    stiffness matrix (multiply by k_cal to get N/m).
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('capsid', pdb_file)

    coords_angstrom = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    coords_angstrom.append(
                        residue['CA'].get_vector().get_array())

    coords_angstrom = np.array(coords_angstrom)
    coords = coords_angstrom * 1e-10   # Angstrom → metres
    N = len(coords)
    cutoff_sq = (cutoff_angstrom * 1e-10) ** 2

    K = lil_matrix((3 * N, 3 * N))
    for i in range(N):
        for j in range(i + 1, N):
            d = coords[i] - coords[j]
            dist_sq = np.dot(d, d)
            if dist_sq <= cutoff_sq:
                for a in range(3):
                    for b in range(3):
                        val = d[a] * d[b] / dist_sq
                        K[3*i+a, 3*j+b] -= val
                        K[3*j+b, 3*i+a] -= val
                        K[3*i+a, 3*i+b] += val
                        K[3*j+a, 3*j+b] += val
    return coords, K.tocsr(), N

# ============================================================
# SPRING CONSTANT CALIBRATION
# ============================================================
def calibrate_spring_exact(R_m, E_Pa, N_nodes, h_m=2e-9,
                           rho_s=rho_shell, rho_fluid=rho_f):
    """
    Calibrate k_cal (N/m) using exact denominator.
    Sets lowest non-trivial ENM mode = f_pred from Eq. (7).
    """
    f_pred = predict_frequency_exact(R_m, E_Pa, h_m,
                                     rho_s=rho_s, rho_fluid=rho_fluid)
    omega_sq = (2.0 * pi * f_pred) ** 2
    denom = exact_denominator(R_m, h_m, rho_s, rho_fluid)
    m_node = denom * (4.0 / 3.0) * pi * R_m**3 / N_nodes
    k_cal = omega_sq * m_node / (xi_10 ** 2)
    return k_cal, f_pred, m_node

# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":

    # Reference system: CCMV
    R_ccmv = 14e-9      # m
    h_ccmv = 2e-9       # m
    E_ccmv = 1.4e8      # Pa

    f_exact = predict_frequency_exact(R_ccmv, E_ccmv, h_ccmv)
    f_approxC3 = predict_frequency_approx(R_ccmv, E_ccmv, C_val=3.0)
    f_approxC1 = predict_frequency_approx(R_ccmv, E_ccmv, C_val=1.0)

    print("CCMV frequency predictions:")
    print(f"  Exact denominator (3*rho_s*h/R + 3*rho_f) : {f_exact/1e9:.3f} GHz")
    print(f"  Parametrized C=3 (rho_s + 3*rho_f)        : {f_approxC3/1e9:.3f} GHz")
    print(f"  Parametrized C=1 (rho_s + rho_f)          : {f_approxC1/1e9:.3f} GHz")

    # Prestress shifts
    P_ccmv = 1.01e6     # Pa (10 atm)
    P_lambda = 5.07e6   # Pa (50 atm)
    P_phi29 = 6.08e6    # Pa (60 atm)

    print("\nFractional prestress frequency shifts:")
    print(f"  CCMV (P=10 atm, E=140 MPa) : {prestress_shift(P_ccmv, E_ccmv)*100:.3f}%")
    print(f"  lambda (P=50 atm, E=100 MPa): {prestress_shift(P_lambda, 1e8)*100:.3f}%")
    print(f"  phi29  (P=60 atm, E=100 MPa): {prestress_shift(P_phi29, 1e8)*100:.3f}%")

    # ENM from PDB (if file supplied)
    if len(sys.argv) > 1:
        pdb_file = sys.argv[1]
        print(f"\nBuilding ENM from {pdb_file} ...")
        coords, K, N = build_enm(pdb_file, cutoff_angstrom=10.0)
        print(f"Number of CA nodes: {N}")

        k_cal, f_pred, m_node = calibrate_spring_exact(
            R_ccmv, E_ccmv, N, h_ccmv)

        print(f"Calibrated spring constant k_cal = {k_cal:.4e} N/m")
        print(f"Effective mass per node m_node = {m_node:.4e} kg")
        print(f"Target frequency (exact denom.) = {f_pred/1e9:.3f} GHz")

        K_scaled = K * k_cal
        D_scaled = K_scaled / m_node

        n_modes = min(30, 3 * N - 6)
        eigenvalues, _ = eigsh(D_scaled, k=n_modes + 6,
                               which='LM', sigma=0)
        eigenvalues = np.sort(np.real(eigenvalues))
        non_rigid = eigenvalues[6:]
        frequencies = np.sqrt(np.abs(non_rigid)) / (2.0 * pi)

        print("\nLowest 10 non-trivial ENM mode frequencies:")
        for i, fq in enumerate(frequencies[:10]):
            print(f"  Mode {i+1:2d}: {fq/1e9:.3f} GHz")

        # Supplementary Figure S1
        plt.figure(figsize=(8, 5))
        plt.plot(range(1, len(frequencies) + 1),
                 frequencies / 1e9, 'o-', markersize=6, color='steelblue')
        plt.xlabel('Mode index', fontsize=12)
        plt.ylabel('Frequency (GHz)', fontsize=12)
        plt.title('CCMV ENM Vibrational Spectrum (PDB: 1CWP)\n'
                  'Calibrated with exact thin-shell denominator', fontsize=11)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig('Figure_S1_ENM_CCMV.pdf', dpi=300)
        plt.savefig('Figure_S1_ENM_CCMV.png', dpi=300)
        print("\nSupplementary Figure S1 saved.")
    else:
        print("\nNo PDB file provided. Run with: python enm_capsid.py 1cwp.pdb")