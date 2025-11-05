from libraries.CrystalStructure import CrystalStructure
import numpy as np

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
                 epsilon: float = 0.345):
        self.crystal = crystal  # oggetto CrystalStructure
        self.sigma   = sigma
        self.epsilon = epsilon

    def compute_potential(self):
            """
            Calcola il potenziale di Lennard-Jones per ogni atomo.
            Il potenziale dipende dai parametri A, B e, indirettamente, dalla distanza di
            taglio R_C, usata per trovare i vicini.
            Restituisce il potenziale totale del sistema.
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
        Versione vettoriale di compute_potential che usa self.crystal.distance_matrix e
        self.crystal.neighbour_matrix. Restituisce il potenziale totale.
        """
        # assicurati che le matrici dei vicini siano presenti
        mask = getattr(self.crystal, "neighbour_matrix", None)
        dist = getattr(self.crystal, "distance_matrix", None)
        if mask is None or dist is None:
            # prova a calcolare i vicini se disponibile la versione numpy
            if hasattr(self.crystal, "find_neighbours_numpy"):
                self.crystal.find_neighbours_numpy(self.crystal.R_C)
                mask = self.crystal.neighbour_matrix
                dist = self.crystal.distance_matrix
            else:
                raise RuntimeError("distance_matrix o neighbour_matrix mancanti; esegui prima find_neighbours_numpy()")

        # estrai distanze solo per le coppie considerate vicine
        r = dist[mask]  # vettore delle distanze r_ij per cui mask è True
        if r.size == 0:
            pot = 0.0
        else:
            s = self.sigma
            e = self.epsilon
            # calcolo vettoriale del LJ
            phi = 4.0 * e * ( (s / r)**12 - (s / r)**6 )
            pot = 0.5 * np.sum(phi)   # fattore 1/2 per evitare doppio conto

        self.crystal.potential = float(pot)
        return self.crystal.potential


    def compute_forces(self):
        """
        Versione che utilizza list.
        Calcola le forze sugli atomi dovute al potenziale di Lennard-Jones (F = -∇V).
        Restituisce una lista di vettori forza per ogni atomo.
        """
        
        def _addendo_derivata_lennard_jones(q_i, q_k, r_ik):
            return 1/(r_ik**8) * ( (2*self.sigma**6)/(r_ik**6) - 1 ) * (q_i - q_k)
            
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
        Versione vettoriale di compute_forces_matrix che usa 
        self.crystal.distance_matrix e self.crystal.neighbour_matrix.
        Restituisce una matrice (N,3) con le forze su ciascun atomo.
        """
        pos = np.asarray(self.crystal.positions, dtype=float)  # (N,3)
        n = pos.shape[0]

        mask = getattr(self.crystal, "neighbour_matrix", None)
        dist = getattr(self.crystal, "distance_matrix", None)

        if mask is None or dist is None:
            if hasattr(self.crystal, "find_neighbours_numpy"):
                self.crystal.find_neighbours_numpy()
                mask = self.crystal.neighbour_matrix
                dist = self.crystal.distance_matrix
            else:
                # fallback: costruisci distanze da pos e usa R_C
                diffs_tmp = pos[:, None, :] - pos[None, :, :]
                dist = np.linalg.norm(diffs_tmp, axis=2)
                mask = (dist <= self.crystal.R_C) & (~np.eye(n, dtype=bool))

        # differenze vettoriali e distanze
        diffs = pos[:, None, :] - pos[None, :, :]   # (N,N,3)
        r = np.asarray(dist, dtype=float)           # (N,N)

        # coefficiente del termine
        coeff = 24.0 * self.epsilon * (self.sigma ** 6)

        # calcola lo scalare per le coppie considerate (evita divisioni per zero)
        scalar = np.zeros_like(r)
        idx = mask
        if np.any(idx):
            r_idx = r[idx]
            inv_r8 = 1.0 / (r_idx ** 8)
            term = (2.0 * (self.sigma ** 6) / (r_idx ** 6)) - 1.0
            scalar[idx] = coeff * inv_r8 * term

        # moltiplica per il vettore differenza e somma sulle colonne j per ottenere la forza su i
        forces = np.sum(scalar[..., None] * diffs, axis=1)  # (N,3)

        return forces
    
    
    def compute_forces_matrix(self):
        """
        Versione che utilizza matrici numpy.
        Calcola le forze sugli atomi dovute al potenziale di Lennard-Jones (F = -∇V).
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