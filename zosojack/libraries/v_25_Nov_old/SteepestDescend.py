import numpy as np
from pathlib import Path
from libraries.CrystalStructure import CrystalStructure
from libraries.CrystalPotential import CrystalPotential
from libraries.PolynomialJunction import PolynomialJunction


class SteepestDescend:
    
    
    def __init__ (self, crystal: CrystalStructure):
        self.crystal = crystal
        
        
    def _output_state(self, filename, step, E_pot, F_max):
        """Scrive lo stato della simulazione (come in CrystalDynamics)."""

        with open(filename, "w") as f:
            f.write(f"{step} \t {E_pot} \t {F_max}")
        
    def _output_positions(self, foldername, step):
            """Salva le posizioni finali (stesso formato dell'originale)."""
            
            filename = foldername / (f"fcc100a{self.crystal.N_atoms}_1.txt")
            np.savetxt(filename, self.crystal.positions, comments="")
            
            
    def minimize_energy(self, max_steps=1000, F_tol=1e-5, C_steep=0.01, pol_junction=False, output=False):
        
        # Predisposizione output
        if output:
            out_dir = Path(
                f"output/min_positions/max_s{max_steps}~C{C_steep}~F_t{F_tol}~Ag~{self.crystal.N_atoms}"
            )
            out_dir.mkdir(parents=True, exist_ok=True)
            state_file = out_dir / "energy.txt"
                        
            n_print = 10
        
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
                forces = potenziale.compute_forces_numpy()
            else:
                potenziale = CrystalPotential(self.crystal)
                forces = potenziale.compute_forces_numpy()
            potential_energy = potenziale.compute_potential_numpy()
            
            max_force = np.max(np.linalg.norm(forces, axis=1))
            if max_force < F_tol:
                print(f"Converged in {step} steps.")
                break
            
            # aggiorno posizioni
            self.crystal.positions += C_steep * forces 
            
            if step == max_steps - 1:
                print("Maximum number of steps reached without convergence.")
                
            
            if output and step % n_print == 0:
                self._output_state(
                    state_file,
                    step,
                    potential_energy,
                    max_force,
                )
        if output:
            self._output_positions(out_dir)