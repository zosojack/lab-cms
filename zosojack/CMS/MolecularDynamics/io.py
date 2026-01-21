# io.py
import numpy as np
from pathlib import Path

np.linalg.norm

from CMS.MolecularDynamics.CrystalStructure import CrystalStructure

class AtomTracker:
    """
    AtomTracker
    ===========
    Classe che traccia le posizioni di uno specifico atomo durante la simulazione.
    
    Attributes
    ----------
    index : int
        Indice dell'atomo da tracciare.
    output_file : str
        Percorso del file di output dove salvare le posizioni.
    pcb_option : str
        Modalità per le condizioni al contorno: 'periodic' oppure 'unbounded' (default).
        
    Methods
    -------
    record_position(step: int, crystal: CrystalStructure) -> None
        Registra la posizione dell'atomo al passo specificato.
    """
    
    def __init__(self, index: int, output_file: str, pcb_option: str = 'unbounded') -> None:
        """Inizializza il tracciatore di atomi con l'indice dell'atomo da tracciare."""

        self.index = index
        self.output_path = Path(output_file) # converte in Path
        self.set_pcb_option(pcb_option)
        
        # Crea la cartella se non esiste
        parent_folder = self.output_path.parent
        parent_folder.mkdir(parents=True, exist_ok=True)
        
        with open(self.output_path, 'w') as f:
            f.write(f"TRAIETTORIA ATOMO {self.index} - {self.pcb_option.upper()}\n")
            f.write(f"step \t x \t y \t z\n")
            
    def set_pcb_option(self, option: str) -> None:
        """
        Imposta l'opzione per il trattamento delle condizioni al contorno periodiche.

        Parameters
        ----------
        option : str
            'periodic' per condizioni al contorno periodiche, 'unbounded' altrimenti.
        """
        if option not in ['periodic', 'unbounded']:
            raise ValueError("L'opzione deve essere 'periodic' o 'unbounded'.")
        self.pcb_option = option

    def record_position(self, step: int, crystal: CrystalStructure) -> None:
        """
        Registra le posizioni degli atomi tracciati nel file di output.

        Parameters
        ----------
        step : int
            Passo temporale corrente della simulazione.
        crystal : CrystalStructure
            Struttura cristallina con posizioni e dimensioni della cella per le PBC.
        """
        pos = crystal.positions[self.index].copy()
        
        if self.pcb_option == 'periodic':
            # Applica le condizioni al contorno periodiche
            pos -= crystal.pcb * np.floor(crystal.positions[self.index] / crystal.pcb)
            
        with open(self.output_path, 'a') as f:
            f.write(f"{step} \t {pos[0]} \t {pos[1]} \t {pos[2]}\n")

class XYZwriter():
    """
    XYZwriter
    ==========
    Classe che registra le traiettorie atomiche in formato .xyz.
    
    Attributes
    ----------
    output_folder : str
        Cartella di output dove salvare le traiettorie.
    dt : float
        Intervallo di tempo tra i frame registrati.
    dump_interval : int
        Passi tra i frame salvati (default: 200).
        
    Methods
    -------
    write_frame(step: int, positions: np.ndarray) -> None
        Registra le posizioni di tutti gli atomi nel file di output.
    """
    
    def __init__(self, output_folder: str, dt: float, dump_interval: int = 200) -> None:
        """ Inizializza l'oggetto che registra le traiettorie degli atomi. """
        self.output_folder = Path(output_folder) # converte in Path
        self.dt = dt
        self.dump_interval = dump_interval
        
        # Crea la cartella se non esiste
        self.output_folder.mkdir(parents=True, exist_ok=True)

    def write_frame(self, step: int, positions: np.ndarray) -> None:
        """
        Registra le posizioni di tutti gli atomi nel file di output.
        L'azione è eseguita soltanto ogni dump_interval steps.

        Parameters
        ----------
        positions : np.ndarray
            Posizioni attuali di tutti gli atomi.
        step : int
            Passo temporale corrente della simulazione.
        """
        if step % self.dump_interval == 0:
            file_index = step // self.dump_interval
            N_atoms = positions.shape[0]
            
            with open(self.output_folder / f"frame_{file_index}.xyz", 'w') as f:
                pass
            
            header = f"{N_atoms}\n time={step * self.dt}"
            np.savetxt(self.output_folder / f"frame_{file_index}.xyz", positions, header=header, comments="")