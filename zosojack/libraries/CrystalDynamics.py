# CrystalDynamics.py
from __future__ import annotations

import warnings

import numpy as np
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, Optional

from libraries.CrystalStructure import CrystalStructure
from libraries.CrystalPotential import CrystalPotential
from libraries.PolynomialJunction import PolynomialJunction
from libraries.io import AtomTracker, XYZwriter

# Costanti fisiche 
k_B = 1 / 11603  # eV/K
N_A = 1.66053906660E-27  # numero di Avogadro

class CrystalDynamics:
	"""
	CrystalDynamics
	===============
	Classe per eseguire la dinamica molecolare classica su un sistema cristallino.
	"""

	def __init__(
		self,
		crystal: CrystalStructure,
		atomic_mass: float = 108,
		dt: float = 1e-15,
		temp_ini: float = 20,
		atom_tracker: Optional[AtomTracker | list[AtomTracker]] = None,
		xyz_writer: Optional[XYZwriter] = None,
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
  
		# Strumenti di output
		self._atom_tracker = atom_tracker
		self._xyz_writer = xyz_writer

	# ---------------------------------------------------------------
	# Properties - non necessarie ma faccio pratica
	# ---------------------------------------------------------------
 
	@property
	def xyz_writer(self):
		return self._xyz_writer

	@xyz_writer.setter
	def xyz_writer(self, writer: XYZwriter) -> None:
		if writer.dump_interval < 10:
			warnings.warn("⚠️ Attenzione: un dump_interval troppo piccolo può rallentare la simulazione e saturare il disco.", UserWarning)
   
		self._xyz_writer = writer 
   	
	@property
	def atom_tracker(self):
		return self._atom_tracker

	@atom_tracker.setter
	def atom_tracker(self, tracker: AtomTracker | list[AtomTracker]) -> None:
		self._atom_tracker = tracker if isinstance(tracker, list) else [tracker]
        
	# ---------------------------------------------------------------
	# Utility private
	# ---------------------------------------------------------------
 
	def set_seed(self, myseed) -> None:
		"""Imposta il seed del generatore pseudo-casuale."""
		np.random.seed(myseed)

	def _random_velocities(self) -> None:
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

	def _update_positions(self) -> None:
		"""Aggiorna le posizioni con Velocity-Verlet (stessa formula)."""

		acc = self.old_force / self.atomic_mass
		displacement = self.velocities * self.dt + 0.5 * acc * (self.dt**2)
		self.crystal.positions += displacement

	def _update_neighbours_distances(self) -> None:
		"""
		Aggiorna le distanze tra gli atomi.
		Se le posizioni non sono cambiate più di un valore soglia, sono ricalcolate soltanto
		le distanze tra i vicini già noti con only_neighbours_distance().
		Se le posizioni sono cambiate troppo, il flag neighbours_computed è posto False e
		viene rieseguito find_neighbours().
		"""	
    	# per aumentare la severità del controllo, uso 0.4*(R_V - R_C) invece di 0.5
		reference_displacement = self.crystal.positions - self.crystal.reference_positions
		if np.max(np.linalg.norm(reference_displacement, axis=1)) > 0.4 * (self.crystal.R_V - self.crystal.R_C):
			self.crystal.neighbours_computed = False
			return self.crystal.only_neighbours_distance()
		else:	
			return self.crystal.find_neighbours()
    
	def _update_forces(self, poly7=None) -> None:
		"""Calcola le forze usando la versione numpy del potenziale."""

		self.new_force = CrystalPotential(self.crystal, poly7=poly7).compute_forces()

	def _update_velocities(self, rescale_velocity: bool = False) -> None:
		"""
  		Aggiorna le velocità.
  		Se rescale_velocity è True, ricalcola la velocità per mantenere la temperatura costante.
  		"""

		acc_old = self.old_force / self.atomic_mass
		acc_new = self.new_force / self.atomic_mass

		self.velocities = self.velocities + 0.5 * (acc_old + acc_new) * self.dt
		if rescale_velocity:
			# Re-scale velocità per mantenere la temperatura costante
			E_K = 0.5 * self.atomic_mass * float(np.sum(self.velocities**2))
			T_now = (2 / 3) * (E_K / (self.crystal.N_atoms * k_B))
			self.velocities *= np.sqrt(self.temp_ini / T_now)
   
		self.kinetic_E = 0.5 * self.atomic_mass * float(np.sum(self.velocities**2))
  
		self.old_force = self.new_force

	def _temperature(self) -> float:
		"""Restituisce la temperatura attuale del sistema."""

		if self.kinetic_E is not None:
			return (2 / 3) * (self.kinetic_E / (self.crystal.N_atoms * k_B))
		else:
			return self.temp_ini

	# NOTE: deprecated, usare CrystalDynamicsResult
	def _output_state(self, filename, step, E_tot, E_pot, E_kin, temp) -> None:
		"""Scrive lo stato della simulazione (come in CrystalDynamics)."""

		with open(filename, "w") as f:
			f.write(f"{step * self.dt} \t {E_tot} \t {E_pot} \t {E_kin} \t {temp}\n")

	# NOTE: deprecated, usare XYZwriter
	def _output_positions(self, foldername, step, n_steps) -> None:
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
                    rescale_velocity: bool = False,
                    debug: bool = False,
                    track_last=False, n_print=None, output=False #NOTE: deprecated arguments
                ) -> CrystalDynamicsResult:
     
		"""Esegue la dinamica molecolare per `n_steps` step."""			

		if getattr(self.crystal, "neighbour_matrix", None) is None:
			print(
				f"⚠️ Vicini non calcolati in precedenza. " + \
        		f"Calcolo con R_C={self.crystal.R_C}, R_P={self.crystal.R_P} e R_V={self.crystal.R_V}."
			)
			self.crystal.find_neighbours()
   
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
			self.old_force = CrystalPotential(self.crystal, poly7=poly7).compute_forces()
	

		# LOOP TERMALIZZAZIONE ------------------------------
		steps_th = int(t_th / self.dt)
		for step in range(steps_th):
			self._update_positions()
			self._update_neighbours_distances()
			self._update_forces(poly7=poly7)
			self._update_velocities(rescale_velocity=rescale_velocity)
		# FINE LOOP TERMALIZZAZIONE -------------------------
	
		# Vettori per salvare i dati
		meta_energies_dict = {
			'total': np.zeros(n_steps - steps_th),
			'kinetic': np.zeros(n_steps - steps_th),
			'potential': np.zeros(n_steps - steps_th),
		}
		meta_T = np.zeros(n_steps - steps_th)
		# NOTE: kinetic_E è sempre aggiornata in _update_velocities()
		# rende tutto un po' confuso 
  		# perché è l'unica energia ad essere un membro invece che una variabile locale
		# LOOP DINAMICA -------------------------------------  
		for step in range(int(n_steps-steps_th)):
			self._update_positions()
			self._update_neighbours_distances()
			self._update_forces(poly7=poly7)
			self._update_velocities()

			potential_energy_now = CrystalPotential(
				self.crystal,
				poly7=poly7,
			).compute_potential()
			temp = self._temperature()

			E_tot_now = potential_energy_now + self.kinetic_E
			meta_energies_dict['total'][step] = E_tot_now
			meta_energies_dict['kinetic'][step] = self.kinetic_E
			meta_energies_dict['potential'][step] = potential_energy_now
			meta_T[step] = temp

			if debug:
				print(
					f"step {step+1}/{n_steps}: E_tot={E_tot_now:.3f} eV, "
					f"V={potential_energy_now:.3f} eV, "
					f"K={self.kinetic_E:.3f} eV, T={temp:.1f} K"
				)
     
			# Tracking adatom specifico
			if self.atom_tracker is not None:
				for a_t in self.atom_tracker:
					a_t.record_position(step, self.crystal)
			# Tracking posizioni
			if self.xyz_writer is not None:
				self.xyz_writer.write_frame(step, self.crystal.positions)
		# FINE LOOP DINAMICA ------------------------------------
  
		# Risultati finali
		result = CrystalDynamicsResult(
			time_step = self.dt,
			num_steps = n_steps,
			steps_th = steps_th,
			energies = meta_energies_dict,
			temperatures = meta_T,
			trajectory_folder_path = self.xyz_writer.output_folder if self.xyz_writer is not None else None,  # Opzionale
			adatom_file_path = [a_t.output_path for a_t in self.atom_tracker] if self.atom_tracker is not None else None, # Opzionale
		)
   
		return result

@dataclass
class CrystalDynamicsResult:
	"""Classe per gestire i risultati di una simulazione di dinamica molecolare."""
 
	# --- 1. CAMPI OBBLIGATORI (vanno passati quando crei l'oggetto) ---
	time_step: float
	num_steps: int
	steps_th: int
	energies: Dict[str, np.ndarray]
	temperatures: np.ndarray

	# --- 2. CAMPI OPZIONALI (hanno un valore di default) ---
	trajectory_folder_path: Optional[str] = None
	adatom_file_path: Optional[str] = None

	# Per liste/dict vuoti di default, si usa field(default_factory=...)
	metadata: Dict = field(default_factory=dict)

	# --- 3. CAMPI CALCOLATI (non li passi, li calcola lui) ---
	# Usiamo init=False per dire "non chiederlo nel costruttore"
	mean_temp:  float = field(init=False)
	final_temp: float = field(init=False)
	mean_E_tot: float = field(init=False)
	std_E_tot:  float = field(init=False)

	def __post_init__(self):
		"""
		Questo metodo viene eseguito automaticamente SUBITO DOPO l'__init__.
		Usato per calcolare le statistiche derivate (medie e dev st).
		"""
		# Calcolo automatico delle temperature
		self.mean_temp = np.mean(self.temperatures)
		self.std_temp = np.std(self.temperatures)
		self.final_temp = self.temperatures[-1]
		
		# Calcolo automatico energia (se presente nel dizionario)
		if 'total' in self.energies:
			self.mean_E_tot = np.mean(self.energies['total'])
			self.std_E_tot = np.std(self.energies['total'])
		else:
			self.mean_E_tot = np.nan
			self.std_E_tot = np.nan
   

	# Metodo per stampare info utili
	def summary(self) -> None:
		print (f"Simulation Result:\n"
				f" - Duration: {self.time_step*self.num_steps:.2f} ps\n"
				f" - Mean Temp: {self.mean_temp:.2f} ± {self.std_temp:.2f} K\n"
				f" - Mean Energy: {self.mean_E_tot:.2f} ± {self.std_E_tot:.2f} eV\n")