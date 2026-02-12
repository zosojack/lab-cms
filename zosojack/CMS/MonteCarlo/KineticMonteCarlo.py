# KineticMonteCarlo.py
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional
from CMS.MonteCarlo.KMCNumbaSubroutines import (
    pbc_corr,
    count_NN,
    find_move,
    _update_NN_deposition,
    _update_NN_diffusion,
)
from CMS.MonteCarlo.ioKMC import XYZwriter

import numpy as np

# Costanti fisiche 
k_B = 1 / 11603  # eV/K

class KineticMonteCarlo:
    """
    KineticMonteCarlo
    =================
    Classe per eseguire simulazioni kinetic Monte Carlo di crescita su reticolo.
    
    Parameters
    ----------
    L : tuple[int, int]
        Dimensioni del reticolo (Lx, Ly).
    flux : float
        Flusso di deposizione in ML/s.
    T : float
        Temperatura in Kelvin.
    nu : float, optional
        Frequenza di tentativi di diffusione (default è 1.E+13 Hz).
    J0 : float, optional
        Energia di interazione del primo vicino (default è 4 * J1 in eV).
    J1 : float, optional
        Energia di interazione del secondo vicino (default è -0.345 eV).
    seed : int, optional
        Seme per il generatore di numeri casuali (default è 123413432).
        
    Attributes
    ----------
    height : np.ndarray
        Matrice che rappresenta l'altezza della superficie in ogni sito del reticolo.
    first_neigh : np.ndarray
        Matrice che contiene il numero di primi vicini per ciascun atomo sulla superficie.
    k_diff : np.ndarray
        Matrice dei rate di diffusione per ciascun sito del reticolo.
    k_diff_sum : float
        Somma totale dei rate di diffusione.
        
    Methods
    -------
    run(end_time: float) -> KineticMonteCarloResult
        Esegue la simulazione fino al tempo specificato e restituisce i risultati
        raggruppati in un oggetto KineticMonteCarloResult.
    """
    
    def __init__(
        self,
        L: tuple[int, int], # (Lx, Ly)
        flux: float,
        T: float,
        nu: float = 1.E+13, # diffusion rate [Hz]
        J0: float = 4*(-0.345), # default: J0 = 4 * J1 [eV]
        J1: float = -0.345, # 2nd neighbor interaction energy [eV]
        xyz_writer: Optional[XYZwriter] = None,
        seed: int = 123413432,
    ):
    
        self.Lx, self.Ly = (L[0],L[1])
        self.flux = flux
        self.T = T if T != 0 else 1E-10  # evita divisioni per zero
        self.nu = nu
        self.J0 = J0
        self.J1 = J1
        self.first_neigh: Optional[np.ndarray] = None
        self.k_diff: Optional[np.ndarray] = None
        self.k_diff_sum: Optional[float] = None
        self.k_diff_row_sums: Optional[np.ndarray] = None
        # strumenti di output
        self.xyz_writer = xyz_writer
        
        np.random.seed(seed)
        
        # inizializzo la matrice delle altezze
        self.height = np.zeros((self.Lx, self.Ly), dtype=np.float64)
        
        # HACK: un atomo può avere 0, 1, 2, 3, 4 primi vicini
        # i inizializzo una volta sola i 5 possibili rate di diffusione
        self.rates_lookup = np.zeros(5, dtype=np.float64)
        inv_kT = 1.0 / (k_B * self.T)
        for n in range(5):
            Eb = -(self.J0 + n * self.J1)
            self.rates_lookup[n] = 4 * self.nu * np.exp(-Eb * inv_kT)

        # inizializzo i rate di diffusione
        self._initialize_diffusion_rates()

        
    def _initialize_diffusion_rates (self) -> None:
        ''' Inizializza i primi vicini, i rate di diffusione e la loro somma '''
        # se non sono stati inizializzati i primi vicini, vanno calcolati
        if self.first_neigh is None:
            self.first_neigh = count_NN(self.height)
        
        # se non sono stati inizializzati i rate di diffusione, vanno calcolati
        if self.k_diff is None:
            # barriera di energia
            Eb = -(self.J0 + self.first_neigh * self.J1)
            # rate di diffusione 
            self.k_diff = 4 * self.nu * np.exp(-Eb / (k_B*self.T))
        
        # se non è stata inizializzata la somma dei rate di diffusione, va calcolata
        if self.k_diff_sum is None:
            self.k_diff_sum = np.sum(self.k_diff)
            
        # se non è stata inizializzata la somma dei rate di diffusione per riga, va calcolata
        if self.k_diff_row_sums is None:
            self.k_diff_row_sums = np.sum(self.k_diff, axis=1)
    
    def _deposition_event (self) -> None:
        # estrae casualmente la posizione in cui depositare un atomo
        x_site = np.random.randint(0, self.Lx)
        y_site = np.random.randint(0, self.Ly)
        # aggiunge un atomo sulla superficie in quella posizone
        self.height[x_site, y_site] += 1.0
        # aggiorna i primi vicini e i rate di diffusione
        self.first_neigh, self.k_diff, self.k_diff_sum, self.k_diff_row_sums = \
            _update_NN_deposition(self.first_neigh,
                                  self.k_diff,
                                  self.k_diff_sum,
                                  self.k_diff_row_sums,
                                  self.height,
                                  (x_site, y_site),
                                  self.rates_lookup,
        )
        
    def _diffusion_event (self,
                          rho: float,
                          k_diff: np.ndarray
    ) -> None:
        # sceglie l'atomo da spostare
        x0, y0 = find_move(rho=rho, 
                           k_diff=k_diff, 
                           current_k_diff_row_sums=self.k_diff_row_sums)
        # sceglie in che direzione spostarlo
        dir = np.random.rand()
        if   dir < 0.25: # sx
            x_d, y_d = pbc_corr(x0-1, self.Lx), y0
        elif dir < 0.50: # dx
            x_d, y_d = pbc_corr(x0+1, self.Lx), y0
        elif dir < 0.75: # up
            x_d, y_d = x0, pbc_corr(y0+1, self.Ly)
        else:          #down
            x_d, y_d = x0, pbc_corr(y0-1, self.Ly)
            
        # decrementa di 1 l'altezza nella posizione di partenza x_0, y_0
        self.height[x0, y0] -= 1.0    
        # incrementa di 1 l'altezza nella posizione di arrivo x_d, y_d
        self.height[x_d, y_d] += 1
        # aggiorna i primi vicini e i rate di diffusione
        self.first_neigh, self.k_diff, self.k_diff_sum, self.k_diff_row_sums = \
            _update_NN_diffusion(self.first_neigh,
                                 self.k_diff,
                                 self.k_diff_sum,
                                 self.k_diff_row_sums,
                                 self.height,
                                 (x0, y0),
                                 (x_d, y_d),
                                 self.rates_lookup,
        )
        
    def _control_array_size (
        self,
        array: np.ndarray,
        current_step: int,
    ) -> np.ndarray:
        if current_step >= len(array):
            # raddoppia la dimensione dell'array
            new_size = len(array) * 2
            new_array = np.empty(new_size, dtype=array.dtype)
            new_array[:len(array)] = array
            return new_array
        else:
            return array
        
    def run (
        self,
        end_time: float,
    ) -> KineticMonteCarloResult:
        # inizializzione variabili
        time = 0.0
        k_depo = self.flux * self.Lx * self.Ly  # tasso di deposizione: prendo il totale poi metto in un sito a caso
        
        # HACK: per evitare append, dichiaro un array con una stima degli eventi di deposizione,
        # poi raddoppio la dimensione se non dovesse bastare
        # Stima accurata: numero atteso di eventi di deposizione
        expected_events = int(self.flux * self.Lx * self.Ly * end_time)
        # inizializzo contenitori risultati
        dt_list = np.empty(expected_events, dtype=np.float64)
        rms_roughness_list = np.empty(expected_events, dtype=np.float64)
        n_diffusion_events = 0
        n_deposition_events = 0
        
        # loop di dinamica KMC
        while time < end_time:
            # calcolo k_tot
            k_tot = k_depo + np.sum(self.k_diff)
        
            # tempo dello step
            tau = - np.log(np.random.rand()) / k_tot
            
            # selezione dell'evento da eseguire
            rho = k_tot * np.random.rand()
            if rho < k_depo:
                self._deposition_event()
                # NOTE: tau e roughness sono registrati solo per eventi di deposizione
                dt_list[n_deposition_events] = time
                rms_roughness_list[n_deposition_events] = np.std(self.height)
                n_deposition_events += 1
                # Tracking deposizioni
                if self.xyz_writer is not None:
                    if time >= self.xyz_writer.start_time: # per ritardare il primo frame se necessario
                        self.xyz_writer.write_frame(time, n_deposition_events, self.height)
            else:
                rho -= k_depo
                self._diffusion_event(rho, self.k_diff)
                n_diffusion_events += 1
            
            # controllo la dimensione degli array
            dt_list = self._control_array_size(dt_list, n_deposition_events)
            rms_roughness_list = self._control_array_size(rms_roughness_list, n_deposition_events)
            # aggiorno il tempo
            time += tau
            
        return KineticMonteCarloResult(
            n_deposition_events=n_deposition_events,
            n_diffusion_events=n_diffusion_events,
            # Slice per prendere solo gli eventi registrati prima di calcolare diff
            dt_list=np.diff(dt_list[:n_deposition_events]), 
            rms_roughness_list=rms_roughness_list[:n_deposition_events],
            k_depo=k_depo,
        )
            
@dataclass
class KineticMonteCarloResult:
    """
    KineticMonteCarloResult
    =======================
    Classe per gestire i risultati di una simulazione di kinetic Monte Carlo.
    """
    # campi forniti in input
    n_deposition_events: int
    n_diffusion_events: int
    dt_list: np.ndarray
    rms_roughness_list: np.ndarray
    k_depo: float
    # campi calcolati in post-init
    scarto_lunghezza_array:  int = field(init=False)
    
    def __post_init__(self):
        self.scarto_lunghezza_array = len(self.dt_list) - self.n_deposition_events
        # taglio gli array alla lunghezza effettiva
        self.dt_list = self.dt_list[:self.n_deposition_events]
        self.rms_roughness_list = self.rms_roughness_list[:self.n_deposition_events]