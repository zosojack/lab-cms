# CrystalDynamics.py
from __future__ import annotations

import warnings

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, Optional
from numba import njit

from CMS.MolecularDynamics.CrystalStructure import CrystalStructure
from CMS.MolecularDynamics.CrystalPotential import CrystalPotential
from CMS.MolecularDynamics.PolynomialJunction import PolynomialJunction
from CMS.MolecularDynamics.io import AtomTracker, XYZwriter

# Costanti fisiche 
k_B = 1 / 11603  # eV/K
N_A = 1.66053906660E-27  # numero di Avogadro

# funzioni jit
@njit(cache=True)
def _remove_rotation(positions: np.ndarray, 
                     velocities: np.ndarray, 
                     atomic_mass: float) -> np.ndarray:
    
    N = positions.shape[0]
    
    # Calcolo del Centro di Massa #
    r_cdm = np.zeros(3, dtype=np.float64)
    for i in range(N):
        r_cdm[0] += positions[i, 0]
        r_cdm[1] += positions[i, 1]
        r_cdm[2] += positions[i, 2]
    r_cdm /= N

    # Inizializzazione accumulatori #
    L_tot = np.zeros(3, dtype=np.float64)
    I = np.zeros((3, 3), dtype=np.float64)
    
    # Loop unico su tutti gli atomi per riempire L e I #
    # Usiamo variabili temporanee scalari per velocità
    for i in range(N):
        # Braccio (distanza dal cdm)
        rx = positions[i, 0] - r_cdm[0]
        ry = positions[i, 1] - r_cdm[1]
        rz = positions[i, 2] - r_cdm[2]
        
        vx = velocities[i, 0]
        vy = velocities[i, 1]
        vz = velocities[i, 2]
        
        # Momento Angolare: L += r x (m*v)
        # Nota: moltiplico per la massa alla fine per risparmiare N operazioni
        L_tot[0] += ry * vz - rz * vy
        L_tot[1] += rz * vx - rx * vz
        L_tot[2] += rx * vy - ry * vx
        
        # Tensore d'Inerzia (I) #
        # Elementi diagonali
        I[0, 0] += ry**2 + rz**2
        I[1, 1] += rx**2 + rz**2
        I[2, 2] += rx**2 + ry**2
        
        # Elementi fuori diagonale (negativi)
        I[0, 1] -= rx * ry
        I[0, 2] -= rx * rz
        I[1, 2] -= ry * rz

    # Applico la massa (scalare) a tutto il tensore e al vettore L
    L_tot *= atomic_mass
    I *= atomic_mass
    
    # Simmetrizzo il tensore d'inerzia
    I[1, 0] = I[0, 1]
    I[2, 0] = I[0, 2]
    I[2, 1] = I[1, 2]
    
    # Calcolo omega #
    # Invece di inv(I) @ L, risolviamo il sistema lineare I * omega = L
    # È più stabile numericamente e piace di più a Numba
    # Gestiamo il caso patologico (determinante nullo) con un try/except non supportato bene in njit puro,
    # quindi assumiamo che il sistema 3D non sia degenere (non lineare).
    omega = np.linalg.solve(I, L_tot)
    
    # Correggo le velocità #
    # v_new = v_old - (omega x r)
    # Modifico l'array velocities in-place o ne creo uno nuovo
    v_corrected = np.zeros_like(velocities)
    
    for i in range(N):
        rx = positions[i, 0] - r_cdm[0]
        ry = positions[i, 1] - r_cdm[1]
        rz = positions[i, 2] - r_cdm[2]
        
        # Prodotto vettoriale omega x r
        wx_r_x = omega[1] * rz - omega[2] * ry
        wx_r_y = omega[2] * rx - omega[0] * rz
        wx_r_z = omega[0] * ry - omega[1] * rx
        
        v_corrected[i, 0] = velocities[i, 0] - wx_r_x
        v_corrected[i, 1] = velocities[i, 1] - wx_r_y
        v_corrected[i, 2] = velocities[i, 2] - wx_r_z
        
    return v_corrected
class CrystalDynamics:
	"""
	CrystalDynamics
	===============
	Classe per eseguire la dinamica molecolare classica su un sistema cristallino.
 
	Attributes
	----------
	crystal : CrystalStructure
		Oggetto CrystalStructure che rappresenta la struttura cristallina.
	atomic_mass : float
		Massa atomica in unità di massa atomica (default: 108 u).
	dt : float
		Passo temporale della simulazione in secondi (default: 1e-15 s).
	temp_ini : float
		Temperatura iniziale in Kelvin (default: 20 K).
	atom_tracker : AtomTracker | list[AtomTracker] | None
		Strumento(i) per tracciare la posizione di specifici atomi (opzionale).
	xyz_writer : XYZwriter | None
		Strumento per salvare le posizioni in formato XYZ (opzionale).
  
	Methods
	-------
	run_dynamics(n_steps: int, t_th: float = 0, rescale_velocity: bool = False, debug: bool = False) \
     -> CrystalDynamicsResult
		Esegue la dinamica molecolare per un numero specificato di passi.
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
			warnings.warn("⚠️ Attenzione: un dump_interval troppo piccolo \
       può rallentare la simulazione e saturare il disco.", UserWarning)
   
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

		# correzione deriva traslazionale (sottraggo la media per ciascun asse)
		vel -= vel.mean(axis=0, keepdims=True)
  
		# correzione rotazione
		vel = _remove_rotation(self.crystal.positions, vel, self.atomic_mass)

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
		""" Salva lo stato della simulazione. """

		with open(filename, "w") as f:
			f.write(f"{step * self.dt} \t {E_tot} \t {E_pot} \t {E_kin} \t {temp}\n")

	# NOTE: deprecated, usare XYZwriter
	def _output_positions(self, foldername, step, n_steps) -> None:
		""" Salva le posizioni istantanee (stesso formato dell'originale). """

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
		"""
  		Esegue la dinamica molecolare per `n_steps` step.
    	
		Parameters
		----------
		n_steps : int
			Numero totale di passi della simulazione.
		t_th : float, optional
			Tempo di termalizzazione in secondi (default: 0 s).
		rescale_velocity : bool, optional
			Se True, ricalcola la velocità per mantenere la temperatura costante (default: False).
		debug : bool, optional
			Se True, stampa informazioni di debug ad ogni passo (default: False).
		track_last : bool, optional
			Se True, traccia l'ultimo atomo (default: False).
		n_print : int | None, optional
			Numero di passi tra le stampe di debug (default: None).
		output : bool, optional
			Se True, salva lo stato della simulazione ad ogni passo (default: False).

		Returns
		-------
		CrystalDynamicsResult
			Oggetto contenente i risultati della simulazione.
     	"""			

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
		steps_th = int(t_th / self.dt) if t_th > 0 else 0
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
				self.xyz_writer.write_frame(step, self.crystal)
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
	"""
 	CrystalDynamicsResult
	=====================
 	Classe per gestire i risultati di una simulazione di dinamica molecolare.

	Attributes
	----------
	time_step : float
		Il passo temporale della simulazione in secondi.
	num_steps : int
		Il numero totale di passi della simulazione.
	steps_th : int
		Il numero di passi di termalizzazione.
	energies : Dict[str, np.ndarray]
		Un dizionario contenente array di energie ('total', 'kinetic', 'potential').
	temperatures : np.ndarray
		Un array contenente le temperature ad ogni passo.
	trajectory_folder_path : Optional[str]
		Il percorso della cartella contenente i file di traiettoria XYZ (se salvati).
	adatom_file_path : Optional[str]
		Una lista di percorsi dei file di tracking degli adatom (se salvati).
	mean_temp : float
		La temperatura media calcolata dopo la simulazione.
	final_temp : float
		La temperatura finale della simulazione.
	mean_E_tot : float
		L'energia totale media calcolata dopo la simulazione.
	std_E_tot : float
		La deviazione standard dell'energia totale calcolata dopo la simulazione.
  
	Methods
	-------
	summary() -> None
		Stampa un riepilogo dei risultati della simulazione.
	"""
 
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
   
	def summary(self) -> None:
		""" Stampa un riepilogo dei risultati della simulazione. """
		print (f"Simulation Result:\n"
				f" - Duration: {self.time_step*self.num_steps:.2f} ps\n"
				f" - Mean Temp: {self.mean_temp:.2f} ± {self.std_temp:.2f} K\n"
				f" - Mean Energy: {self.mean_E_tot:.2f} ± {self.std_E_tot:.2f} eV\n")