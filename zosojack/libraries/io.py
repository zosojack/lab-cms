# AtomTracker.py
'''
------------------------
Classe che traccia gli atomi di un sistema cristallino durante una simulazione di dinamica.
------------------------
'''

from pathlib import Path
import numpy as np

# FIXME: import circolare
'''from libraries.CrystalDynamics import CrystalDynamics

def _output_name (dynamics: CrystalDynamics, n_steps: int) -> None:

    out_dir = Path(
        f"output/dynamics/steps{dynamics.n_steps}~dt{dynamics.dt}~T{dynamics.temp_ini}~Ag~{dynamics.crystal.N_atoms}"
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    state_file = out_dir / "energy.txt"
    # TODO: sistemare tracking adatoms in modo intelligente
    track_file = out_dir / "adatom_track.txt"
'''

class AtomTracker:
    
    def __init__(self, index: int, output_file: str) -> None:
        """
        Inizializza il tracciatore di atomi con l'indice dell'atomo da tracciare.
        
        Parametri:
        - index: indice dell'atomo da tracciare
        - output_file: percorso del file di output dove salvare le posizioni
        """
        self.index = index
        self.output_path = Path(output_file) # converte in Path
        
        # Crea la cartella se non esiste
        parent_folder = self.output_path.parent
        parent_folder.mkdir(parents=True, exist_ok=True)
        
        with open(self.output_path, 'w') as f:
            f.write(f"TRAIETTORIA ATOMO {self.index}\n")
            f.write(f"step \t x \t y \t z\n")

    def record_position(self, step: int, positions: np.ndarray) -> None:
        """
        Registra le posizioni degli atomi tracciati nel file di output.
        
        Parametri:
        - step: passo temporale corrente della simulazione
        - positions: array delle posizioni attuali di tutti gli atomi
        """
        with open(self.output_path, 'a') as f:
            pos = positions[self.index]
            f.write(f"{step} \t {pos[0]} \t {pos[1]} \t {pos[2]}\n")

class XYZwriter():
    
    def __init__(self, output_folder: str, dt: float, dump_interval: int) -> None:
        """
        Inizializza l'oggetto che registra le traiettorie degli atomi.
        
        Parametri:
        - output_folder: percorso della cartella di output dove salvare le traiettorie
        - dt: intervallo di tempo tra i frame registrati
        - dump_interval: intervallo di passi tra i frame registrati
        """
        self.output_folder = Path(output_folder) # converte in Path
        self.dt = dt
        self.dump_interval = dump_interval
        
        # Crea la cartella se non esiste
        self.output_folder.mkdir(parents=True, exist_ok=True)

    def write_frame(self, step: int, positions: np.ndarray) -> None:
        
        """
        Registra le posizioni di tutti gli atomi nel file di output.
        L'azione è eseguita soltanto ogni dump_interval steps.
        
        Parametri:
        - positions: array delle posizioni attuali di tutti gli atomi
        - step: passo temporale corrente della simulazione
        """
        if step % self.dump_interval == 0:
            file_index = step // self.dump_interval
            N_atoms = positions.shape[0]
            
            with open(self.output_folder / f"frame_{file_index}.xyz", 'w') as f:
                pass
            
            header = f"{N_atoms}\n time={step * self.dt}"
            np.savetxt(self.output_folder / f"frame_{file_index}.xyz", positions, header=header, comments="")