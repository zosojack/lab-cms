# SteepestDescent.py
import numpy as np
from pathlib import Path

from libraries.CrystalStructure import CrystalStructure
from libraries.CrystalPotential import CrystalPotential
from libraries.PolynomialJunction import PolynomialJunction

class SteepestDescend:
    """
    SteepestDescent
    ===============
    Classe che implementa la minimizzazione dell'energia tramite il metodo del gradiente discendente.
    """
    
    def __init__ (self, crystal: CrystalStructure):
        self.crystal = crystal
        
            
    def minimize_energy(self, max_steps=1000, F_tol=1e-5, C_steep=0.01, pol_junction=False):
        
        if getattr(self.crystal, "neighbour_matrix", None) is None:
            print(
                f"⚠️ Vicini non calcolati in precedenza. " + \
                f"Calcolo con R_C={self.crystal.R_C} e R_P={self.crystal.R_P}."
            )
            self.crystal.find_neighbours()
        
        # Array forza massima e energia
        max_forces_list = []
        potential_energies_list = []
        
        for step in range (max_steps):
            
            # calcolo forze
            if pol_junction:
                poly7 = PolynomialJunction(
                    R_C=self.crystal.R_C,
                    R_P=self.crystal.R_P,  # esempio di punto di giunzione
                    epsilon=0.345,  # eV
                    sigma=2.644,  # Å
                )
                potenziale = CrystalPotential(self.crystal, poly7=poly7)
                forces = potenziale.compute_forces()
            else:
                potenziale = CrystalPotential(self.crystal)
                forces = potenziale.compute_forces()
            potential_energy = potenziale.compute_potential()
                    
            max_force = np.max(np.linalg.norm(forces, axis=1))
            
            potential_energies_list.append(potential_energy)
            max_forces_list.append(max_force)
            
            if max_force < F_tol:
                print(f"Converged in {step} steps.")
                break
            
            # aggiorno posizioni
            self.crystal.positions += C_steep * forces 
            # aggiorno vicini e distanze
            self.crystal.find_neighbours()
            
            if step == max_steps - 1:
                print("Maximum number of steps reached without convergence.")
            
            
        return np.array(potential_energies_list), np.array(max_forces_list)