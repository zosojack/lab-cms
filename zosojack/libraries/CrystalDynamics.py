"""CrystalDynamics.py
========================

Versione di `CrystalDynamics` che sfrutta le versioni vettoriali
(`find_neighbours_numpy`, `compute_potential_numpy`, `compute_forces_numpy`).
Il comportamento esterno rimane identico all'implementazione
originale basata su liste.
"""

import numpy as np
from pathlib import Path

from libraries.CrystalStructure import CrystalStructure
from libraries.CrystalPotential import CrystalPotential
from libraries.PolynomialJunction import PolynomialJunction


# Costanti fisiche 
k_B = 1 / 11603  # eV/K
N_A = 1.66053906660E-27  # numero di Avogadro

class CrystalDynamics:
	"""Versione numba della dinamica molecolare classica."""

	def __init__(
		self,
		crystal: CrystalStructure,
		atomic_mass: float = 108,
		dt: float = 1e-15,
		temp_ini: float = 20,
	):
		# Oggetti
		self.crystal = crystal  # oggetto CrystalStructure

		# Parametri della simulazione (stessa formula del codice originale)
		self.atomic_mass = atomic_mass * 1.66e-27 / 16  # massa in eV s^2/Å^2
		self.dt = dt  # passo temporale
		self.temp_ini = temp_ini  # temperatura iniziale

		# Variabili di stato
		self.old_force = None
		self.new_force = None
		self.velocities = None
		self.temp_final = None
		self.kinetic_E = None
		self.potential_E = None

	# ---------------------------------------------------------------
	# Utility private
	# ---------------------------------------------------------------
	def set_seed(self, myseed):
		"""Imposta il seed del generatore pseudo-casuale."""

		np.random.seed(myseed)

	def _random_velocities(self):
		"""Inizializza le velocità casuali (stesso algoritmo dell'originale)."""

		N_atoms = self.crystal.N_atoms
		C = np.sqrt((3 * k_B * self.temp_ini) / self.atomic_mass)
		vel = np.random.uniform(-C, C, size=(N_atoms, 3))

		# correzione drift (sottraggo la media per ciascun asse)
		vel -= vel.mean(axis=0, keepdims=True)

		# energia cinetica prima del rescaling (identica a prima)
		E_K = 0.5 * self.atomic_mass * float(np.sum(vel**2))
		T_prime = (2 / 3) * (E_K / (N_atoms * k_B))

		self.velocities = vel * np.sqrt(self.temp_ini / T_prime)
		self.kinetic_E = E_K

	def _update_positions(self):
		"""Aggiorna le posizioni con Velocity-Verlet (stessa formula)."""

		acc = self.old_force / self.atomic_mass
		self.crystal.positions += (
			self.velocities * self.dt + 0.5 * acc * (self.dt**2)
		)

	def _update_neighbours(self):
		"""Aggiorna i vicini usando la versione numpy del metodo."""

		self.crystal.find_neighbours_numba()

	def _update_forces(self, poly7=None):
		"""Calcola le forze usando la versione numpy del potenziale."""

		self.new_force = CrystalPotential(self.crystal, poly7=poly7).compute_forces_numba()

	def _update_velocities(self):
		"""Aggiorna le velocità (formula invariata)."""

		acc_old = self.old_force / self.atomic_mass
		acc_new = self.new_force / self.atomic_mass

		self.velocities = self.velocities + 0.5 * (acc_old + acc_new) * self.dt
		self.kinetic_E = 0.5 * self.atomic_mass * float(np.sum(self.velocities**2))
  
		self.old_force = self.new_force

	def _temperature(self):
		"""Restituisce la temperatura corrente del sistema."""

		if self.kinetic_E is not None:
			return (2 / 3) * (self.kinetic_E / (self.crystal.N_atoms * k_B))
		else:
			return self.temp_ini

	def _output_state(self, filename, step, E_tot, E_pot, E_kin, temp):
		"""Scrive lo stato della simulazione (come in CrystalDynamics)."""

		with open(filename, "w") as f:
			f.write(f"{step * self.dt} \t {E_tot} \t {E_pot} \t {E_kin} \t {temp}")

	def _output_positions(self, foldername, step, n_steps):
		"""Salva le posizioni istantanee (stesso formato dell'originale)."""

		width = len(str(n_steps))
		step_str = f"{step:0{width}d}"

		filename = foldername / (f"t={self.dt}~" + step_str + f"_{n_steps}.xyz")
		header = f"{self.crystal.N_atoms}\n time={step * self.dt}"
		np.savetxt(filename, self.crystal.positions, header=header, comments="")

	# ---------------------------------------------------------------
	# Ciclo di dinamica
	# ---------------------------------------------------------------
	def run_dynamics(self, 
                  	n_steps: float,
                    t_th:float = 0, 
                    output: bool = False,
                    debug: bool = False):
     
		"""Esegue la dinamica molecolare per `n_steps` step."""

		if output:
			out_dir = Path(
				f"output/dynamics/steps{n_steps}~dt{self.dt}~T{self.temp_ini}~Ag~{self.crystal.N_atoms}"
			)
			out_dir.mkdir(parents=True, exist_ok=True)
			state_file = out_dir / "energy.txt"
			n_print = 10

		if getattr(self.crystal, "neighbour_matrix", None) is None:
			print(
				f"⚠️ Vicini non calcolati in precedenza. " + \
        		f"Calcolo con R_C={self.crystal.R_C} e R_P={self.crystal.R_P}."
                
			)
			self.crystal.find_neighbours_numba()
   
		if np.isfinite(self.crystal.R_P):
			poly7 = PolynomialJunction(
				R_C=self.crystal.R_C,
				R_P=self.crystal.R_P,  # punto di giunzione
				epsilon=0.345,  # eV
				sigma=2.644,  # Å
			)
		else:
			poly7 = None

		if self.velocities is None:
			self._random_velocities()

		if self.old_force is None:
			self.old_force = CrystalPotential(self.crystal, poly7=poly7).compute_forces_numba()
		
		meta_E_tot = []
		meta_E_K = []
		meta_T = []

		# LOOP TERMALIZZAZIONE
		steps_th = int(t_th / self.dt)
		for step in range(steps_th):
			self._update_positions()
			self._update_neighbours()
			self._update_forces(poly7=poly7)
			self._update_velocities()
   
		for step in range(steps_th, n_steps):
			self._update_positions()
			self._update_neighbours()
			self._update_forces(poly7=poly7)
			self._update_velocities()

			potential_energy = CrystalPotential(
				self.crystal,
				poly7=poly7,
			).compute_potential_numba()
			temp = self._temperature()

			E_tot_now = potential_energy + self.kinetic_E
			meta_E_tot.append(E_tot_now)
			meta_E_K.append(self.kinetic_E)
			meta_T.append(temp)

			if debug:
				print(
					f"step {step+1}/{n_steps}: E_tot={E_tot_now:.3f} eV, "
					f"V={potential_energy:.3f} eV, "
					f"K={self.kinetic_E:.3f} eV, T={temp:.1f} K"
				)

			if output and step % n_print == 0:
				self._output_state(
					state_file,
					step,
					E_tot_now,
					potential_energy,
					meta_E_K,
					meta_T,
				)
				self._output_positions(out_dir, step, n_steps)

		return meta_E_tot, meta_E_K, meta_T