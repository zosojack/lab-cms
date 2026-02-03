# stampare TUTTE LE ENERGIE
# non necessario per ora, si può tenere in RAM

# stampare la matrice delle altezze SOLO quando l'energia è minore di quella precedente
# cioè quando trovo un nuovo minimo

# ioKMC.py
import numpy as np
from pathlib import Path
from numba import njit

@njit
def _convert_pile_to_xyz(height: np.ndarray) -> np.ndarray:
    ''' 
    Converte una matrice di altezze (pile di atomi) in una matrice di posizioni (N_atoms x 3).
    Ogni riga della matrice risultante rappresenta le coordinate (x, y, z) di un atomo.
    '''
    Lx, Ly = height.shape
    N_atoms = int(np.sum(height)) + Lx * Ly  # include gli atomi del substrato a z=0
    positions = np.zeros((N_atoms, 3))
    
    index = 0
    for x in range(Lx):
        for y in range(Ly):
            # Aggiungo atomo del substrato a z=0
            positions[index, 0] = x
            positions[index, 1] = y
            positions[index, 2] = 0
            index += 1
            
            # Aggiungo atomi della pila da z=1 a h
            h = int(height[x, y])
            for z in range(1, h + 1):
                positions[index, 0] = x
                positions[index, 1] = y
                positions[index, 2] = z
                index += 1
                
    return positions

class XYZwriter:
    """ 
    XYZwriter
    =========
    
    Classe per registrare le traiettorie atomiche in formato XYZ. 
    """
    
    def __init__(self, 
                 output_folder: str,
                 max_files: int = 500) -> None:
        """ Inizializza l'oggetto che registra le traiettorie degli atomi. """
        self.output_folder = Path(output_folder) # converte in Path
        self.max_files = max_files
        self.current_file_count = 0
        
        # Crea la cartella se non esiste
        self.output_folder.mkdir(parents=True, exist_ok=True)
        
    
    def write_frame(self, current_accepted_step: int, height: np.ndarray) -> None:
        '''
        Costruisce e registra le posizioni degli atomi depositati sulla superficie.
        L'azione è eseguita un massimo di max_files volte.

        Parameters
        ----------
        current_accepted_step : int
            Numero di volte che è stata accettata una mossa della simulazione.
        height : np.ndarray
            Altezze delle pile di atomi sulla superficie.
        '''
    
        # da pile di atomi, devo passare a posizioni (N_atomsxN_atoms -> N_atomsx3)
        positions = _convert_pile_to_xyz(height)
        
        if self.current_file_count < self.max_files:
            
            N_atoms = positions.shape[0]
            
            with open(self.output_folder / f"frame_{self.current_file_count}.xyz", 'w') as f:
                pass
            
            header = f"{N_atoms}\n time={current_accepted_step:.1f}"
            np.savetxt(self.output_folder / f"frame_{self.current_file_count}.xyz", positions, header=header, comments="")
            
            self.current_file_count += 1
            
        else:
            return