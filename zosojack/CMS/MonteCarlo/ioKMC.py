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
                 max_deposition_steps: int = 1_000, 
                 max_files: int = 500,
                 start_time: float = 0.0) -> None:
        """ Inizializza l'oggetto che registra le traiettorie degli atomi. """
        self.output_folder = Path(output_folder) # converte in Path
        self.max_deposition_steps = max_deposition_steps
        self.max_files = max_files
        self.current_file_count = 0
        self.start_time = start_time # to delay the first frame if needed
        self.dump_interval = max(self.max_deposition_steps // max_files, 1)

        # HACK
        self.start_deposition_step = None
        
        # Crea la cartella se non esiste
        self.output_folder.mkdir(parents=True, exist_ok=True)
        
    
    def write_frame(self, time: float, deposition_step: int, height: np.ndarray) -> None:
        '''
        Costruisce e registra le posizioni degli atomi depositati sulla superficie.
        L'azione è eseguita un massimo di max_files volte.

        Parameters
        ----------
        time : float
            Tempo corrente della simulazione.
        deposition_step : int
            Passo temporale corrente della simulazione.
        height : np.ndarray
            Altezze delle pile di atomi sulla superficie.
        '''
        # HACK
        self.start_deposition_step = deposition_step if self.start_deposition_step is None else self.start_deposition_step
        
        current_writing_step = deposition_step - self.start_deposition_step
        
        # early exit
        if current_writing_step > self.max_deposition_steps:
            return
        
        # da pile di atomi, devo passare a posizioni (N_atomsxN_atoms -> N_atomsx3)
        positions = _convert_pile_to_xyz(height)
        
        if current_writing_step % self.dump_interval == 0:
            
            file_index = int(current_writing_step // self.dump_interval)
            N_atoms = positions.shape[0]
            
            with open(self.output_folder / f"frame_{file_index}.xyz", 'w') as f:
                pass
            
            header = f"{N_atoms}\n time={time:.4f}"
            np.savetxt(self.output_folder / f"frame_{file_index}.xyz", positions, header=header, comments="")
            
            self.current_file_count += 1

            
    """   
    def write_frame(self, time: float, current_dep_step: int, height: np.ndarray) -> None:
        '''
            Appende un frame al file XYZ se le condizioni sono soddisfatte.
        '''
        
        # 1. CONTROLLO "EARLY EXIT": Se abbiamo superato il 10% del layer, smettiamo di scrivere.
        #    La simulazione continua, ma noi non intasiamo l'hard disk.
        if current_dep_step > self.max_deposition_steps:
            return

        # 2. CONTROLLO INTERVALLO: Salviamo solo ogni tot passi
        if current_dep_step % self.dump_interval == 0:
            
            # Convertiamo le altezze in posizioni (N_atoms x 3)
            # Nota: se height è vuoto (step 0), positions sarà vuoto. Gestiamolo.
            positions = _convert_pile_to_xyz(height)
            N_atoms = positions.shape[0]
            
            # Se non ci sono atomi (step 0), possiamo saltare o scrivere un dummy header
            if N_atoms == 0:
                return

            # Costruiamo l'header richiesto da Ovito/XYZ standard
            # Linea 1: Numero atomi
            # Linea 2: Commento (usiamo il tempo)
            header = f"{N_atoms}\nProperties=pos:R:3 Time={time:.4f}"
            
            # Apriamo in modalità APPEND ('a')
            with open(self.filepath, 'a') as f:
                # Scriviamo l'header manualmente
                f.write(header + "\n")
                # Usiamo numpy per scrivere la matrice di coordinate in modo efficiente
                np.savetxt(f, positions, fmt='%.0f')"""