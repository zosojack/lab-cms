# CrystalPotential.py
import numpy as np
from numba import njit

from libraries.CrystalStructure import CrystalStructure
from libraries.PolynomialJunction import PolynomialJunction

@njit(cache=True, fastmath=True)
def _potential_kernel(dist_matrix, neighbour_mask, second_mask, sigma, epsilon, coeffs):
    pot = 0.0
    N = dist_matrix.shape[0]
    for i in range(N):
        for j in range(i + 1, N):  # conti solo coppie (i<j) per evitare 1/2
            r = dist_matrix[i, j]
            if not np.isfinite(r) or r == 0.0:
                continue

            if neighbour_mask[i, j]:
                inv_r = 1.0 / r
                inv_r6 = (sigma * inv_r) ** 6
                inv_r12 = inv_r6 * inv_r6
                pot += 4.0 * epsilon * (inv_r12 - inv_r6)
            elif second_mask[i, j] and coeffs is not None:
                A,B,C,D,E,F,G,H = coeffs
                pot += A + B*r + C*r**2 + D*r**3 + E*r**4 + F*r**5 + G*r**6 + H*r**7
    return pot
    
@njit(cache=True)
def _forces_kernel(positions, dist_matrix, neighbour_mask, second_mask, sigma, epsilon, coeffs):
    N = positions.shape[0]
    forces = np.zeros((N, 3), dtype=np.float64)
    s6 = sigma**6

    for i in range(N):
        xi0, yi0, zi0 = positions[i, 0], positions[i, 1], positions[i, 2]
        for j in range(N):
            if i == j:
                continue
            r = dist_matrix[i, j]
            if not np.isfinite(r) or r == 0.0:
                continue

            dx = xi0 - positions[j, 0]
            dy = yi0 - positions[j, 1]
            dz = zi0 - positions[j, 2]

            if neighbour_mask[i, j]:
                inv_r = 1.0 / r
                inv_r2 = inv_r * inv_r
                inv_r6 = inv_r2 * inv_r2 * inv_r2
                term = (2.0 * s6 * inv_r6) - 1.0
                inv_r8 = inv_r6 * inv_r2
                Fscalar = 24.0 * epsilon * s6 * term * inv_r8
                forces[i, 0] += Fscalar * dx
                forces[i, 1] += Fscalar * dy
                forces[i, 2] += Fscalar * dz
            elif second_mask[i, j] and coeffs is not None:
                # derivata polinomio: B + 2Cr + 3Dr^2 + ... + 7Hr^6
                B,C,D,E,F,G,H = coeffs[1], coeffs[2], coeffs[3], coeffs[4], coeffs[5], coeffs[6], coeffs[7]
                r2 = r * r
                r3 = r2 * r
                r4 = r3 * r
                r5 = r4 * r
                r6 = r5 * r
                dVdr = B + 2*C*r + 3*D*r2 + 4*E*r3 + 5*F*r4 + 6*G*r5 + 7*H*r6
                Fscalar = -dVdr / r
                forces[i, 0] += Fscalar * dx
                forces[i, 1] += Fscalar * dy
                forces[i, 2] += Fscalar * dz

    return forces
    
class CrystalPotential:
    """
    Classe per calcolare il potenziale e le forze in una struttura cristallina.
    Attributi:
    - crystal: oggetto CrystalStructure
    - sigma, epsilon: parametri del potenziale di Lennard-Jones,
                      di default sono quelli per l'argento
    """
    def __init__(self, 
                 crystal: CrystalStructure, 
                 sigma: float = 2.644, 
                 epsilon: float = 0.345,
                 poly7: PolynomialJunction = None):
        self.crystal = crystal  # oggetto CrystalStructure
        self.sigma   = sigma
        self.epsilon = epsilon
        self.poly7 = poly7  # polinomio di settimo grado per la giunzione polinomiale
    # ---------------------------------------------------------------
        
    def compute_potential_numba(self):
        mask_first = self.crystal.neighbour_matrix
        mask_second = self.crystal.second_neighbour_matrix
        dist = self.crystal.distance_matrix

        coeffs = None
        if self.poly7 is not None:
            coeffs = np.array([self.poly7.A, self.poly7.B, self.poly7.C, self.poly7.D,
                            self.poly7.E, self.poly7.F, self.poly7.G, self.poly7.H], dtype=np.float64)

        pot = _potential_kernel(dist, mask_first, mask_second, self.sigma, self.epsilon, coeffs)
        self.crystal.potential = float(pot)
        return self.crystal.potential

    def compute_forces_numba(self):
        pos = np.asarray(self.crystal.positions, dtype=np.float64)
        dist = self.crystal.distance_matrix
        m1 = self.crystal.neighbour_matrix
        m2 = self.crystal.second_neighbour_matrix

        coeffs = None
        
        if self.poly7 is None and np.isinf(self.crystal.R_P) is False:
                raise ValueError("PolynomialJunction non definito per le seconde vicine.")
        elif self.poly7 is not None:
            coeffs = np.array([self.poly7.A, self.poly7.B, self.poly7.C, self.poly7.D, self.poly7.E, self.poly7.F, self.poly7.G, self.poly7.H], dtype=np.float64)

        return _forces_kernel(pos, dist, m1, m2, self.sigma, self.epsilon, coeffs)