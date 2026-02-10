# MetropolisMonteCarlo.py
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from CMS.MonteCarlo.ioMMC import XYZwriter
from CMS.MonteCarlo.MMCNumbaSubroutines import _compute_system_energy, _update_system_energy

# Costanti fisiche 
k_B = 1 / 11603  # eV/K

class MetropolisMonteCarlo:
    """
    MetropolisMonteCarlo
    ====================
    Classe per eseguire simulazioni Monte Carlo di tipo Metropolis su reticolo. 
    Può individuare la configurazione di equilibrio di un sistema di atomi su un reticolo 2D, 
    non prevederne le traiettorie. Non contiene variabili dipendenti dal tempo.
    
    Parameters
    ----------
    L : tuple[int, int]
        Dimensioni del reticolo (Lx, Ly).
    N_atoms : int
        Numero totale di atomi nel reticolo.
    T : float
        Temperatura in Kelvin.
    J : float
        Energia di interazione tra spin (default è 1.0 eV).
    seed : int, optional
        Seme per il generatore di numeri casuali (default è 123413432).
        
    Attributes
    ----------
    spins : np.ndarray
        Matrice che rappresenta gli spin in ogni sito del reticolo.
        
    Methods
    -------
    run(steps: int) -> MetropolisMonteCarloResult
        Esegue la simulazione per un numero specificato di passi e restituisce i risultati
        raggruppati in un oggetto MetropolisMonteCarloResult.
    """
    
    def __init__(
        self,
        L: tuple[int, int], # (Lx, Ly)
        N_atoms: int,
        T: float,
        J: float =  -0.345,
        seed: int = 123413432,
        xyz_writer: Optional[XYZwriter] = None,
    ):
        
        self.Lx, self.Ly = L[0], L[1]
        self.N_atoms = N_atoms
        self.T = T if T > 0 else 1e-10  # evita divisioni per zero
        self.J = J
        self.seed = seed
        # Inizializza la matrice delle posizioni
        self.height = np.zeros(L, dtype=int)
        # Riempie casualmente la matrice con N_atoms atomi
        for _ in range(N_atoms):
            x = np.random.randint(0, self.Lx)
            y = np.random.randint(0, self.Ly)
            self.height[x, y] += 1
        # Inizializza la lista dei siti occupati # NOTE: prima era matrice numpy
        rows, cols = np.where(self.height > 0)
        self.occupied_sites = list(zip(rows, cols))
        # Inizializza l'energia del sistema
        self.energy = _compute_system_energy(self.height)
        # Inizializza l'energia minima trovata
        self.min_energy = self.energy
        
        # Imposta il generatore di numeri casuali
        np.random.seed(self.seed)
        
        # Impostazioni per output XYZ
        self.xyz_writer = xyz_writer
            
    
    def _select_start_site(self) -> tuple[int, int]:
        ''' Seleziona un sito casuale tra quelli occupati '''
        idx = np.random.randint(0, len(self.occupied_sites))
        return self.occupied_sites[idx]
    
    def _select_end_site(self) -> tuple[int, int]:
        ''' Seleziona un sito casuale tra tutti i siti del reticolo '''
        x = np.random.randint(0, self.Lx)
        y = np.random.randint(0, self.Ly)
        return (x, y)
        
    
            
    def run(self, N_steps: int = 100_000, thermalization_steps: int = 0) -> MetropolisMonteCarloResult:
        
        # array energie
        energies = np.empty(N_steps-thermalization_steps)
        energies[0] = self.energy
        # contatore accettazioni
        accepted_moves = 0

        # loop principale    
        for step in range(1, N_steps):
            # Seleziona il sito di partenza
            start_site = self._select_start_site()
            x0, y0 = start_site
            # Seleziona il sito di arrivo
            end_site = self._select_end_site()
            x1, y1 = end_site
            # CONTROLLA SE LO SPOSTAMENTO È VALIDO #
            # la trial move è effettuata direttamente sulla matrice originale
            self.height[x0, y0] -= 1
            self.height[x1, y1] += 1
            # TODO: AGGIORNARE ENERGIA IN MODO PIÙ EFFICIENTE
            # trial_energy = _compute_system_energy(self.height) # energia trial
            trial_energy = _update_system_energy(
                height=self.height,
                current_energy=self.energy,
                start_site=(x0, y0), # La Z non serve, la legge dentro
                end_site=(x1, y1)
            )
            delta_E = trial_energy - self.energy
            # Se delta_E <= 0 accetta lo spostamento
            if delta_E <= 1E-8: # invece di 0, per evitare problemi numerici
                accepted = True
            # Altrimenti, si estrae un numero casuale e si confronta con un fattore di Boltzmann    
            else:
                r = np.random.rand()
                accepted = r < np.exp(-delta_E / (k_B * self.T))
                
            # Esegui lo spostamento se accettato
            if accepted:
                # la matrice delle altezze è già aggiornata
    
                # aggiorna l'energia del sistema
                self.energy = trial_energy
                # aggiorna la lista dei siti occupati
                if self.height[x0, y0] == 0: # se ora non ce ne sono più, va rimosso lo start site
                    self.occupied_sites.remove((x0, y0))
                if self.height[x1, y1] == 1: # se ora ce n'è proprio uno, va aggiunto l'end site
                    self.occupied_sites.append((x1, y1))
                # se l'energia è diminuita, salvo la configurazione    
                if self.min_energy > self.energy:
                    self.min_energy = self.energy
                    if self.xyz_writer is not None: 
                        self.xyz_writer.write_frame(
                            current_accepted_step=accepted_moves,
                            height=self.height,
                        )
                accepted_moves += 1 if step >= thermalization_steps else 0
            else:
                # se la mossa è rifiutata, va annullata la trial move
                self.height[x0, y0] += 1
                self.height[x1, y1] -= 1

            if step >= thermalization_steps:
                energies[step - thermalization_steps] = self.energy
        
        result = MetropolisMonteCarloResult(
            energies=energies,
            acceptance_ratio=accepted_moves / (N_steps - thermalization_steps),
        )
        
        return result

@dataclass    
class MetropolisMonteCarloResult:
    """
    MetropolisMonteCarloResult
    ==========================
    Classe per memorizzare i risultati di una simulazione Metropolis Monte Carlo.
    
    Parameters
    ----------
    energies : np.ndarray
        Lista delle energie totali ad ogni passo.
    acceptance_ratio : float
        Rapporto di accettazione degli spostamenti proposti.
    """
    
    energies: np.ndarray
    acceptance_ratio: float
    
    min_energy: float = field(init=False)
    avg_energy: float = field(init=False)
    std_energy: float = field(init=False)
        
    def __post_init__(self):
        self.min_energy = np.min(self.energies)
        self.avg_energy = np.mean(self.energies)
        self.std_energy = np.std(self.energies)