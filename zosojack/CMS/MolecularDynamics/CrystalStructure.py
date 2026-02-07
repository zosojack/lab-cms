# CrystalStructure.py
from __future__ import annotations
import psutil
import warnings

import numpy as np
from numba import njit
import copy


# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
# KERNELS DEI METODI DELLA CLASSE CrystalStructure
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

@njit(cache=True)
def _compute_displacement_matrix(positions) -> np.ndarray:
    N = positions.shape[0]
    # N x N x 3: matrice dei vettori spostamento
    disp = np.zeros((N, N, 3), dtype=np.float64)
    
    for i in range(N):
        x_i, y_i, z_i = positions[i, 0], positions[i, 1], positions[i, 2]
        # Ottimizzazione: loop solo sul triangolo superiore
        for j in range(i + 1, N):
            if i == j:
                continue
            dx = x_i - positions[j, 0] 
            dy = y_i - positions[j, 1]
            dz = z_i - positions[j, 2]
            
            # Riempie [i, j]
            disp[i, j, 0] = dx
            disp[i, j, 1] = dy
            disp[i, j, 2] = dz
            
            # Riempie [j, i] con segno opposto (antisimmetria)
            disp[j, i, 0] = -dx
            disp[j, i, 1] = -dy
            disp[j, i, 2] = -dz
            
    return disp

def _compute_distance_matrix(positions) -> np.ndarray:
    disp = _compute_displacement_matrix(positions)
    dist = np.linalg.norm(disp, axis=2)
    return dist

@njit(cache=True)
def _find_neighbour_masks_kernel(positions: np.ndarray,
                                 R_C: float,
                                 R_V: float,
                                 pbc: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    '''
    HACK: poiché rij viene salvato solo per i vicini, i confronti vengono effettuati
    senza calcolarne la radice quadrata, per efficienza. Per farlo, anche R_P, R_C e R_V
    vanno considerati al quadrato.
    '''
    # Pre-calcolo quadrati
    R_C_squared = R_C*R_C
    R_V_squared = R_V*R_V
    # Per i loop
    N = positions.shape[0]
    # Nuove matrici dei vicini e delle distanze
    neighbour_list = np.full((N, N), -1, dtype=np.int32)
    num_neighbours = np.zeros(N, dtype=np.int32)
    dist = np.full((N, N), np.inf) # inf se non vicini, 0 se stessi
    
    
    if pbc is not None:
        Lx, Ly, Lz = pbc[0], pbc[1], pbc[2]
    
    # RICALCOLO TUTTE LE DISTANZE
    for i in range(N):
        x_i, y_i, z_i = positions[i, 0], positions[i, 1], positions[i, 2]
        for j in range(i + 1, N):
            # semplice distanza tra atomi i e j
            dx = x_i - positions[j, 0] # x_i - x_j 
            dy = y_i - positions[j, 1]
            dz = z_i - positions[j, 2]
            # considero anche la periodicità al contorno
            dx -= Lx * np.round(dx / Lx)
            dy -= Ly * np.round(dy / Ly)
            dz -= Lz * np.round(dz / Lz)
            rij_squared = dx*dx + dy*dy + dz*dz
            # se rij è nella cella di Verlet R_V, lo considero
            if rij_squared <= R_V_squared:
                # aggiungo l'indice j alla lista dei vicini di i
                neighbour_list[i, num_neighbours[i]] = j
                # aggiorno il numero di vicini di i
                num_neighbours[i] += 1
                
                # aggiungo l'indice i alla lista dei vicini di j (Simmetria)
                idx_j = num_neighbours[j]
                if idx_j < N:
                    neighbour_list[j, idx_j] = i
                    num_neighbours[j] += 1
                    
                # se rij è anche entro R_C, salvo la distanza
                if rij_squared <= R_C_squared:
                    dist[i, j] = dist[j, i] = np.sqrt(rij_squared)
                
    for k in range(N):
        dist[k, k] = 0.0
        
    return neighbour_list, num_neighbours, dist

@njit(cache=True)
def _only_neighbours_distance_kernel(positions: np.ndarray, 
                                     R_C: float, 
                                     neighbour_list: np.ndarray, 
                                     num_neighbours: np.ndarray, 
                                     pbc: np.ndarray) -> np.ndarray:
    
    '''

    '''
    # Per efficienza, confronto i quadrati delle distanze
    R_C_squared = R_C*R_C
    # Per la periodicità al contorno
    Lx, Ly, Lz = pbc[0], pbc[1], pbc[2]
    # Per i loop
    N = positions.shape[0]
    
    # Nuova matrice delle distanze
    new_dist = np.full((N, N), np.inf) # inf se non vicini, 0 se stessi
    
    # RICALCOLO SOLO LE DISTANZE TRA I VICINI GIÀ NOTI
    # per ciascun atomo i
    for i in range(N):
        x_i, y_i, z_i = positions[i, 0], positions[i, 1], positions[i, 2]
    
        # itero sui suoi vicini j, k-esimi vicini in neighbour_list
        for k in range(num_neighbours[i]):
            j = neighbour_list[i, k]
            
            # Evito doppio calcolo (poiché la lista contiene sia i->j che j->i)
            # Calcolo solo se j > i
            if j <= i:
                continue

            dx = x_i - positions[j, 0] 
            dy = y_i - positions[j, 1]
            dz = z_i - positions[j, 2]
            
            dx -= Lx * np.round(dx / Lx)
            dy -= Ly * np.round(dy / Ly)
            dz -= Lz * np.round(dz / Lz)
            
            rij_squared = dx*dx + dy*dy + dz*dz
            
            # Se la distanza aggiornata è ancora entro R_C, la salvo
            if rij_squared <= R_C_squared:
                d_val = np.sqrt(rij_squared)
                new_dist[i, j] = new_dist[j, i] = d_val
                    
    for k in range(N):
        new_dist[k, k] = 0.0
        
    return new_dist

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
# CLASS CrystalStructure
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
class CrystalStructure:
    """
    CrystalStructure
    ================
    Classe per rappresentare una struttura cristallina.
    
    Attributes
    ----------
    positions : np.ndarray
        Matrice Nx3 delle posizioni atomiche.
    N_atoms : int
        Numero totale di atomi.
    R_C : float
        Distanza di taglio per i primi vicini.
    R_P : float
        Punto di giunzione polinomiale che separa primi e secondi vicini.
    R_V : float
        Raggio della Verlet cage entro cui tenere traccia degli atomi.
    neighbours_computed : bool
        Flag che indica se matrici di vicini e distanze sono aggiornate.
    reference_positions : np.ndarray | None
        Copia delle posizioni quando i vicini sono stati calcolati.
    pbc : np.ndarray
        Dimensioni della cella per le condizioni periodiche.
        
    Properties
    ----------
    vec_x : np.ndarray
        Vettore delle coordinate x degli atomi.
    vec_y : np.ndarray
        Vettore delle coordinate y degli atomi.
    vec_z : np.ndarray
        Vettore delle coordinate z degli atomi.
    displacements_matrix : np.ndarray
        Tensore NxNx3 dei vettori spostamento fra tutti gli atomi.
    crystal_center : np.ndarray
        Coordinate (x, y, z) del centro del volume del cristallo.
    
    Methods
    -------
    from_file(filename) -> CrystalStructure
        Crea un oggetto CrystalStructure leggendo da file.
    empty(n_atoms) -> CrystalStructure
        Costruttore alternativo: crea un oggetto CrystalStructure vuoto.
    copy() -> CrystalStructure
        Restituisce una copia di un oggetto CrystalStructure.
    set_R_C(R_C) -> None
        Imposta la distanza di taglio R_C.
    set_R_P(R_P) -> None
        Imposta il punto di giunzione polinomiale R_P.
    set_R_V(R_V) -> None
        Imposta la grandezza della Verlet cage R_V.
    set_pbc(pbc) -> None
        Imposta la dimensione della cella per la periodicità al contorno.
    add_atom(position) -> None
        Aggiunge un atomo alla struttura cristallina.
    find_neighbours() -> None
        Ricalcola tutte le distanze e trova primi e secondi vicini.
    only_neighbours_distance() -> None
        Ricalcola solo le distanze dai vicini già noti.
    print_neighbours(index=None) -> None
        Stampa gli indici dei primi vicini per ogni atomo o di uno specifico.
    print_second_neighbours(index=None) -> None
        Stampa gli indici dei secondi vicini per ogni atomo o di uno specifico.
    """
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
    # COSTRUTTORI
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
    
    def __init__(self, positions):
        self.positions = positions  # matrice Nx3 delle posizioni
        self.N_atoms = positions.shape[0]  # numero totale di atomi
         
        self.R_C = np.inf  # distanza di taglio usata per trovare i vicini
        self.R_P = np.inf  # punto di giunzione polinomiale
        self.R_V = np.inf  # grandezza della Verlet cage: vicini che non contribuiscono
        self.neighbours_computed = False # flag per calcolo di vicini e distanze
        self.reference_positions = None  # posizioni di riferimento per la Verlet cage
        
        # 16.6416 è la dimensione ideale della cella per Ag in Å
        self.pbc = np.full(3, 1664.16) # *100 così di default non considera la periodicità
        
        # NUOVE STRUTTURE DATI
        # Sostituiscono le matrici booleane N x N
        self.neighbour_list = None     # Matrice (N, MAX_NEIGH) di int
        self.num_neighbours = None     # Array (N) di int (conteggio vicini per ogni atomo)
        self.distance_matrix = None    # Matrice densa distanze per compatibilità
        
    @classmethod
    def from_file(cls, filename) -> CrystalStructure:
        """ Crea un Crystal leggendo da file. """
        data = np.loadtxt(filename)
        return cls(data)
    # Uso:
    # crystal = Crystal.from_file('data.txt')
    
    @classmethod
    def empty(cls, n_atoms) -> CrystalStructure:
        """ Costruttore alternativo: crea Crystal vuoto """
        zeros = np.zeros((n_atoms, 3), dtype=np.float64)
        return cls(zeros)
    
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
    # PROPERTIES
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
    
    @property
    def vec_x(self) -> np.ndarray:
        """ Restituisce il vettore delle coordinate x degli atomi. """
        return self.positions[:, 0]

    @property
    def vec_y(self) -> np.ndarray:
        """ Restituisce il vettore delle coordinate y degli atomi. """
        return self.positions[:, 1]

    @property
    def vec_z(self) -> np.ndarray:
        """ Restituisce il vettore delle coordinate z degli atomi. """
        return self.positions[:, 2]
    
    @property
    def displacements_matrix(self) -> np.ndarray:
        """
        Restituisce il tensore NxNx3 dei vettori spostamento fra tutti gli atomi.
        Shape: (N_atomi, N_atomi, 3)
        """
        N = self.positions.shape[0]
        
        # --- Stima Memoria ---
        bytes_per_float = 8  # float64
        num_elements = N * N * 3
        required_memory_bytes = num_elements * bytes_per_float
        
        # Ottieni la RAM disponibile attuale
        available_ram = psutil.virtual_memory().available
        
        # Conversione per leggibilità (GB)
        req_gb = required_memory_bytes / (1024**3)
        avail_gb = available_ram / (1024**3)

        # Logica di warning (es. se occupiamo più del 20% della RAM libera)
        if required_memory_bytes > (available_ram * 0.2):
            warnings.warn(
                f"\nATTENZIONE: La matrice richiesta occuperà {req_gb:.2f} GB "
                f"su {avail_gb:.2f} GB disponibili. Potrebbe causare swap o crash.",
                ResourceWarning
            )

        # Chiama la funzione numba (definita fuori dalla classe per performance)
        return _compute_displacement_matrix(self.positions)
    
    @property
    def crystal_center(self) -> np.ndarray:
        """
        Restituisce le coordinate (x, y, z) del centro del volume del cristallo.
        """
        
        # passo reticolare 'a' (distanza minima tra primi vicini)
        if getattr(self, "distance_matrix", None) is None:
            dist = _compute_distance_matrix(self.positions)
        else:
            dist = self.distance_matrix
            
        # Maschera gli zeri e trova il minimo globale (scalare)
        passo_reticolare = np.min(np.where(dist < 0.5, np.inf, dist))
        
        min_pos = np.min(self.positions, axis=0)
        max_pos = np.max(self.positions, axis=0)
        PM = (min_pos + max_pos) / 2.0
        
        # punto medio tra questo e quello dopo
        center = PM + (passo_reticolare / 2.0) 
        
        return center
        
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
    # METODI 
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
    
    def copy(self) -> CrystalStructure:
        return copy.deepcopy(self)
    
    def set_R_C(self, R_C) -> None:
        """
        Imposta la distanza di taglio R_C usata per trovare i vicini.

        Parameters
        ----------
        R_C : float
            Distanza di taglio per i primi vicini.
        """
        self.R_C = R_C
        
    def set_R_P(self, R_P) -> None:
        """
        Imposta il punto di giunzione polinomiale R_P; 
        divide primi e 'secondi' vicini.

        Parameters
        ----------
        R_P : float
            Raggio che separa primi e secondi vicini.
        """
        self.R_P = R_P
    
    def set_R_V(self, R_V) -> None:
        """
        Imposta la grandezza della Verlet cage R_V.
        Consigliato: R_V = R_C + 0.5 Å.
        Non può essere minore di R_C, altrimenti solleva un ValueError.

        Parameters
        ----------
        R_V : float
            Raggio della Verlet cage (deve essere >= R_C).
        """
        if R_V < self.R_C:
            raise ValueError(f"R_V ({R_V}) non può essere minore di R_C ({self.R_C}).")
        self.R_V = R_V
        
    def set_pbc(self, pbc: tuple | np.ndarray) -> None:
        """
        Imposta la dimensione della cella per la condizione di periodicità al contorno.
        Deve essere un array-like di 3 elementi.
        Deve essere consistente con R_C.
        Può ricevere np.inf per indicare nessuna periodicità in una direzione.

        Parameters
        ----------
        pbc : array-like
            Dimensioni della cella (Lx, Ly, Lz) in Å; usare np.inf per direzione non periodica.
        """
        # controllo lunghezza
        if len(pbc) != 3:
            raise ValueError("pbc deve essere un array-like di 3 elementi: (Lx, Ly, Lz).")
        # se va bene, converto in np.ndarray di float64
        pbc = np.asarray(pbc, dtype=np.float64)
        # controllo consistenza con R_C
        if self.R_C >= 0.5 * min(pbc):
            raise ValueError("⚠️ Attenzione: R_C deve essere minore della metà della dimensione della cella.")
        
        # per ricevere inf, deve convertirlo in un float molto grande
        # sostituisco prima gli inf
        finite_values = pbc[np.isfinite(pbc)]
        if finite_values.size == 0:
            raise ValueError("Almeno una dimensione della cella deve essere finita.")

        replacement = 100 * np.max(finite_values)
        pbc = np.where(np.isinf(pbc), replacement, pbc)
        
        self.pbc = pbc
    
    # TODO: metodo per aggiungere più atomi contemporaneamente?
    def add_atom(self, position) -> None:
        """
        Aggiunge un atomo alla struttura cristallina.

        Parameters
        ----------
        position : array-like
            Coordinate (x, y, z) dell'atomo da aggiungere.
        """
        position = np.asarray(position, dtype=np.float64)
        if position.shape != (3,):
            raise ValueError("La posizione deve essere un array-like di 3 elementi.")
        self.positions = np.vstack([self.positions, position])
        self.N_atoms += 1
        self.neighbours_computed = False  # i vicini devono essere ricalcolati
    
    def find_neighbours(self) -> None:
        '''
        Ricalcola OGNI distanza e trova primi e secondi vicini per ogni atomo in base a R_C e R_P.
        Gli atomi entro una distanza di Verlet (R_C < r < R_V) sono aggiunti alla matrice dei secondi vicini,
        ma la loro distanza non è salvata in distance_matrix.

        Returns
        -------
        None
        '''
        positions = np.asarray(self.positions, dtype=np.float64)
        neighbour_list, num_neighbours, dist = _find_neighbour_masks_kernel(
            positions=positions,
            R_C=self.R_C, 
            R_V=self.R_V, 
            pbc=self.pbc
        )
        # aggiorno le matrici dei vicini e delle distanze
        self.neighbour_list = neighbour_list
        self.num_neighbours = num_neighbours
        self.distance_matrix = dist
        # eseguito il calcolo, si attiva un flag e si salvano le posizioni di riferimento
        self.neighbours_computed = True
        self.reference_positions = np.copy(self.positions)
        
    def only_neighbours_distance(self) -> None:
        '''
        Ricalcola solamente le distanze dai vicini già noti, senza aggiornarli.
        Se le posizioni sono cambiate troppo, neighbours_computed è posto False ed
        è necessario rieseguire find_neighbours().

        Returns
        -------
        None
        '''
        positions = np.asarray(self.positions, dtype=np.float64)
        dist = _only_neighbours_distance_kernel(
            positions=positions, 
            R_C=self.R_C, 
            neighbour_list=self.neighbour_list, 
            num_neighbours=self.num_neighbours, 
            pbc=self.pbc
        )
        # aggiorno le matrici dei vicini e delle distanze
        self.distance_matrix = dist
        
    # NOTE: deprecated: matrici sostituite da neighbour_list e num_neighbours
    def print_neighbours(self, index=None) -> None:
        """
        Stampa gli indici dei primi vicini per ogni atomo o di uno nello specifico.
        """
        if not self.neighbours_computed:
            print("Prima di chiamare questo metodo, esegui find_neighbours()")
            return None
        
        if index is not None:
            # basta dire quali sono i True sulla riga index di neighbour_matrix:
            print(f"Vicini dell'atomo {index}: {np.where(self.neighbour_matrix[index])[0]}")
            return None
        else:
            print("Indici dei vicini per ogni atomo:")
            for i in range(self.N_atoms):
                print(f"Atomo {i}, n_neigh={np.sum(self.neighbour_matrix[i])}: \
                    {np.where(self.neighbour_matrix[i])[0]}")
    
    # NOTE: deprecated: matrici sostituite da neighbour_list e num_neighbours
    def print_second_neighbours(self, index=None) -> None:
        """
        Stampa gli indici dei secondi vicini per ogni atomo o di uno nello specifico.
        """
        if not self.neighbours_computed:
            print("Prima di chiamare questo metodo, esegui find_neighbours()")
            return None 
        if np.isfinite(self.R_P) is False:
            print("R_P non definito; assegnarlo con CrystalStructure.set_R_P ed eseguire find_neighbours().")
            return None
        if index is not None:
            print(f"Secondi vicini dell'atomo {index}: \
                {np.where(self.second_neighbour_matrix[index])[0]}")
            return None
        else:
            print("Indici dei secondi vicini per ogni atomo:")
            for i in range(self.N_atoms):
                print(f"Atomo {i}, n_2nd_neigh={np.sum(self.second_neighbour_matrix[i])}: \
                    {np.where(self.second_neighbour_matrix[i])[0]}")