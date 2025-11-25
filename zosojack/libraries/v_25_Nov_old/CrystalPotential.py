# CrystalPotential.py

from libraries.CrystalStructure import CrystalStructure
from libraries.PolynomialJunction import PolynomialJunction
import numpy as np
from numba import jit

from numba import njit
    
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
        

    def compute_potential(self):
            """
            Calcola il potenziale di Lennard-Jones per ogni atomo.
            Il potenziale dipende dai parametri A, B e, indirettamente, dalla distanza di
            taglio R_C, usata per trovare i vicini.
                LJ    r < R_P
            V = P7    R_P < r < R_C --> polynomial junction
                0     r > R_C
            Restituisce il potenziale totale del sistema.
            Tipicamente: R_C = 4.5 Å, R_P = 4.2 Å
            """
            def _lennard_jones_ij(r_ij):
                twelve = (self.sigma/r_ij)**12
                six = (self.sigma/r_ij)**6
                return 4*self.epsilon*(twelve - six)
            
            potenziale = 0
            for i, atom_neighbours in enumerate(self.crystal.which_neighbour):
                for j, neighbour_index in enumerate(atom_neighbours):
                    if i != neighbour_index:
                        d_ij = self.crystal.how_distant[i][j]
                        potenziale += _lennard_jones_ij(d_ij)
                
            self.crystal.potential = potenziale / 2
            return self.crystal.potential
        
    def compute_potential_numpy(self):
        """
        Versione vettoriale di compute_potential che usa self.crystal.distance_matrix,
        self.crystal.neighbour_matrix e self.crystal.second_neighbour_matrix.
        Restituisce il potenziale totale.
        """
        # assicurati che le matrici dei vicini siano presenti
        mask_first  = getattr(self.crystal, "neighbour_matrix", None)
        mask_second = getattr(self.crystal, "second_neighbour_matrix", None)
        dist = getattr(self.crystal, "distance_matrix", None)
        if mask_first is None or dist is None or mask_second is None:
            # prova a calcolare i vicini se disponibile la versione numpy
            if hasattr(self.crystal, "find_neighbours_numpy"):
                self.crystal.find_neighbours_numpy(self.crystal.R_C)
                mask_first = self.crystal.neighbour_matrix
                mask_second = self.crystal.second_neighbour_matrix
                dist = self.crystal.distance_matrix
            else:
                raise ValueError("distance_matrix o neighbour_matrix mancanti; esegui prima find_neighbours_numpy()")

        # estrai distanze solo per le coppie considerate PRIME vicine
        r_first = dist[mask_first]  # vettore delle distanze r_ij per cui mask_first è True
        r_second = dist[mask_second]  # vettore delle distanze r_ij per cui mask_second è True
        # calcola il potenziale di Lennard-Jones per le prime vicine
        if r_first.size == 0:
            pot = 0.0
        else:
            s = self.sigma
            e = self.epsilon
            # calcolo vettoriale del LJ
            phi = 4.0 * e * ( (s / r_first)**12 - (s / r_first)**6 )
            pot = 0.5 * np.sum(phi)   # fattore 1/2 per evitare doppio conto
        # calcola il potenziale polinomiale per le seconde vicine
        if r_second.size > 0:
            if self.poly7 is None:
                self.poly7 = PolynomialJunction(epsilon=self.epsilon,
                                                sigma=self.sigma,
                                                R_C=self.crystal.R_C,
                                                R_P=self.crystal.R_P)
            phi_poly = self.poly7.eval(r_second)
            pot += 0.5 * np.sum(phi_poly)  # fattore 1/2 per evitare doppio conto

        self.crystal.potential = float(pot)
        return self.crystal.potential


    def compute_forces(self):
        """
        Versione che utilizza list.
        Calcola le forze sugli atomi dovute al potenziale di Lennard-Jones (F = -∇V).
        Restituisce una lista di vettori forza per ogni atomo.
        # NOTE: NON considera la giunzione polinomiale
        """
        
        def _addendo_derivata_lennard_jones(q_i, q_k, r_ik):
            return 1/(r_ik**8) * ( (2*self.sigma**6)/(r_ik**6) - 1 ) * (q_i - q_k)
        
        def addendo_derivata_polinomio(r_ik, B, C, D, E, F, G, H):
            return -(B + 2*C*r_ik + 3*D*r_ik**2 + 4*E*r_ik**3 + 5*F*r_ik**4 + 6*G*r_ik**5 + 7*H*r_ik**6)
            
        vec_forza = [] # ciascun entrata è un atomo, ciascun atomo è una lista di tre entrate
        for i, atom_neighbours in enumerate(self.crystal.which_neighbour):
            forza_x, forza_y, forza_z = 0, 0, 0
            for j, neighbour_index in enumerate(atom_neighbours):
                d_ij = self.crystal.how_distant[i][j]
                forza_x += _addendo_derivata_lennard_jones(self.crystal.vec_x[i], self.crystal.vec_x[neighbour_index], d_ij)
                forza_y += _addendo_derivata_lennard_jones(self.crystal.vec_y[i], self.crystal.vec_y[neighbour_index], d_ij)
                forza_z += _addendo_derivata_lennard_jones(self.crystal.vec_z[i], self.crystal.vec_z[neighbour_index], d_ij)
            forza_x *= 24 * self.epsilon * self.sigma**6
            forza_y *= 24 * self.epsilon * self.sigma**6
            forza_z *= 24 * self.epsilon * self.sigma**6
            vec_forza.append([forza_x, forza_y, forza_z])
            
        return vec_forza
    
    def compute_forces_numpy(self):
        """
        Usa neighbour_matrix (prime vicine, LJ) e second_neighbour_matrix (giunzione polinomiale).
        Restituisce una matrice (N,3) con le forze su ciascun atomo.
        """
        pos = np.asarray(self.crystal.positions, dtype=float)  # (N,3)
        n = pos.shape[0]

        mask_first  = getattr(self.crystal, "neighbour_matrix", None)
        mask_second = getattr(self.crystal, "second_neighbour_matrix", None)
        dist = getattr(self.crystal, "distance_matrix", None)

        if mask_first is None or mask_second is None or dist is None:
            if hasattr(self.crystal, "find_neighbours_numpy"):
                # usa la tua versione numpy (senza argomenti se così l’hai definita)
                self.crystal.find_neighbours_numpy()
                mask_first  = self.crystal.neighbour_matrix
                mask_second = self.crystal.second_neighbour_matrix
                dist = self.crystal.distance_matrix
            else:
                # fallback minimale
                diffs_tmp = pos[:, None, :] - pos[None, :, :]
                dist = np.linalg.norm(diffs_tmp, axis=2)
                eye = np.eye(n, dtype=np.bool_)
                mask_first = (dist <= self.crystal.R_P) & (~eye)
                mask_second = (dist > self.crystal.R_P) & (dist <= self.crystal.R_C) & (~eye)

        # differenze vettoriali e distanze
        diffs = pos[:, None, :] - pos[None, :, :]   # (N,N,3)
        r = np.asarray(dist, dtype=float)           # (N,N)

        # accumulatore scalare per coppie (i,j); poi verrà moltiplicato per diffs e sommato su j
        scalar = np.zeros_like(r)

        # 1) contributo Lennard-Jones sulle prime vicine
        if np.any(mask_first):
            s6 = self.sigma ** 6
            coeff_LJ = 24.0 * self.epsilon * s6
            r_first = r[mask_first]
            inv_r8 = 1.0 / (r_first ** 8)
            term = (2.0 * s6 / (r_first ** 6)) - 1.0
            scalar[mask_first] += coeff_LJ * inv_r8 * term  # questo è g(r)/r per LJ

        # 2) contributo polinomiale sulle seconde vicine: g(r) = -dV/dr, poi dividi per r
        if np.any(mask_second):
            r_second = r[mask_second]

            # costruisci il polinomio se serve
            if self.poly7 is None and np.isinf(self.crystal.R_P) is False:
                raise ValueError("PolynomialJunction non definito per le seconde vicine.")

            # valuta -dV/dr dal polinomio
            g = self.poly7.eval_derivative(r_second)  # atteso: g(r) = -dV/dr
            
            scalar[mask_second] -= g / r_second  # converte in fattore per (pos_i - pos_j)

        # forza su i: somma su j [ scalar_ij * (pos_i - pos_j) ]
        forces = np.sum(scalar[..., None] * diffs, axis=1)  # (N,3)
        return forces
    
    
    def compute_forces_matrix(self):
        """
        Versione che utilizza matrici numpy.
        Calcola le forze sugli atomi dovute al potenziale di Lennard-Jones (F = -∇V).
        # NOTE: NON considera la componente dovuta alla giunzione polinomiale.
        Restituisce una matrice Nx3 di forze per ogni atomo.
        """
        
        def addendo_derivata_lennard_jones(q_i, q_k, r_ik):
            return 1/(r_ik**8) * ( (2*self.sigma**6)/(r_ik**6) - 1 ) * (q_i - q_k)
            
        mat_forza = np.zeros((self.crystal.N_atoms, 3)) # ciascuna riga è un atomo, ognuna con tre colonne
        for i, atom_neighbours in enumerate(self.crystal.which_neighbour):
            for j, neighbour_index in enumerate(atom_neighbours):
                d_ij = self.crystal.how_distant[i][j]
                mat_forza[i, 0] += addendo_derivata_lennard_jones(self.crystal.vec_x[i], self.crystal.vec_x[neighbour_index], d_ij)
                mat_forza[i, 1] += addendo_derivata_lennard_jones(self.crystal.vec_y[i], self.crystal.vec_y[neighbour_index], d_ij)
                mat_forza[i, 2] += addendo_derivata_lennard_jones(self.crystal.vec_z[i], self.crystal.vec_z[neighbour_index], d_ij)
        mat_forza *= 24 * self.epsilon * self.sigma**6

        return mat_forza
    
    
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
            
    

