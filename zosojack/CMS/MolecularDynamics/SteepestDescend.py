# SteepestDescent.py
import numpy as np
from pathlib import Path

from CMS.MolecularDynamics.CrystalStructure import CrystalStructure
from CMS.MolecularDynamics.CrystalPotential import CrystalPotential
from CMS.MolecularDynamics.PolynomialJunction import PolynomialJunction

class SteepestDescend:
    """
    SteepestDescent
    ===============
    Classe che implementa la minimizzazione dell'energia tramite il metodo del gradiente discendente.
    
    Attributes
    ----------
    crystal : CrystalStructure
        Struttura cristallina da minimizzare.
        
    Methods
    -------
    minimize_energy(max_steps=1000, F_tol=1e-5, C_steep=0.01, pol_junction=False) -> Tuple[np.ndarray, np.ndarray]
        Esegue la minimizzazione dell'energia e restituisce l'energia potenziale e le forze massime ad ogni passo.
    """
    
    def __init__ (self, crystal: CrystalStructure):
        """Inizializza l'oggetto SteepestDescend con la struttura cristallina da minimizzare."""
        self.crystal = crystal
        
            
    def minimize_energy(self, max_steps=1000, F_tol=1e-5, C_steep=0.01, pol_junction=False):
        """
        Esegue la minimizzazione dell'energia e restituisce l'energia potenziale e le forze massime ad ogni passo.

        Parameters
        ----------
        max_steps : int, optional
            Numero massimo di passi per la minimizzazione (default è 1000).
        F_tol : float, optional
            Tolleranza per la forza massima per la convergenza (default è 1e-5).
        C_steep : float, optional
            Coefficiente di discesa per l'aggiornamento delle posizioni (default è 0.01).
        pol_junction : bool, optional
            Se True, utilizza il giunto polinomiale per il potenziale (default è False).
        """
        
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