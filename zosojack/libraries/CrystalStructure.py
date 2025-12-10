# CrystalStructure.py 
import numpy as np
from numba import njit
import copy

@njit(cache=True)
def _find_neighbour_masks_kernel(positions,
                                 R_P,
                                 R_C,
                                 R_V,
                                 pcb):
    '''
    HACK: poiché rij viene salvato solo per i vicini, i confronti vengono effettuati
    senza calcolarne la radice quadrata, per efficienza. Per farlo, anche R_P, R_C e R_V
    vanno considerati al quadrato.
    '''
    # Pre-calcolo quadrati
    R_P_squared = R_P*R_P
    R_C_squared = R_C*R_C
    R_V_squared = R_V*R_V
    # Per i loop
    N = positions.shape[0]
    # Nuove matrici dei vicini e delle distanze
    neighbour = np.zeros((N, N), dtype=np.bool_)
    second = np.zeros((N, N), dtype=np.bool_)
    dist = np.full((N, N), np.inf) # inf se non vicini, 0 se stessi
    
    if pcb is not None:
        Lx, Ly, Lz = pcb[0], pcb[1], pcb[2]
    
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
                # se è maggiore di R_C, lo segno come 2° vicino ma non salvo la distanza
                if rij_squared > R_C_squared:
                    second[i, j] = second[j, i] = True
                # se è minore di R_C, salvo la distanza e controllo se è 1° o 2° vicino
                else:
                    dist[i, j] = dist[j, i] = np.sqrt(rij_squared)
                    # se è entro R_P, lo segno come 1° vicino
                    if rij_squared <= R_P_squared:
                        neighbour[i, j] = neighbour[j, i] = True
                    # se è entro R_C, lo segno come 2° vicino
                    else:
                        second[i, j] = second[j, i] = True
                    
    for k in range(N):
        dist[k, k] = 0.0
    return neighbour, second, dist

@njit(cache=True)
def _only_neighbours_distance_kernel(positions, 
                                     R_P, 
                                     R_C, 
                                     R_V,
                                     neighbour, 
                                     second, 
                                     pcb):
    
    '''
    Docstring:
    '''
    # Per efficienza, confronto i quadrati delle distanze
    R_P_squared = R_P*R_P
    R_C_squared = R_C*R_C
    R_V_squared = R_V*R_V
    # Per la periodicità al contorno
    Lx, Ly, Lz = pcb[0], pcb[1], pcb[2]
    # Per i loop
    N = positions.shape[0]
    # Nuova matrice delle distanze
    new_dist = np.full((N, N), np.inf) # inf se non vicini, 0 se stessi
    
    # RICALCOLO SOLO LE DISTANZE TRA I VICINI GIÀ NOTI
    for i in range(N):
        x_i, y_i, z_i = positions[i, 0], positions[i, 1], positions[i, 2]
        for j in range(i + 1, N):
            # check veloce per vedere se i e j sono vicini
            if not (neighbour[i, j] or second[i, j]):
                continue 
            
            # semplice distanza tra atomi i e j
            dx = x_i - positions[j, 0] # x_i - x_j 
            dy = y_i - positions[j, 1]
            dz = z_i - positions[j, 2]
            # considero anche la periodicità al contorno
            dx -= Lx * np.round(dx / Lx)
            dy -= Ly * np.round(dy / Ly)
            dz -= Lz * np.round(dz / Lz)
            rij_squared = dx*dx + dy*dy + dz*dz
            
            # per la logica di 'vicinanza'
            is_neigh = False
            is_second = False
            d_val = np.inf

            # se rij è ancora nella cella di Verlet R_V, lo considero
            if rij_squared <= R_V_squared:
                # se è maggiore di R_C, lo segno come 2° vicino ma non salvo la distanza
                if rij_squared > R_C_squared:
                    is_second = True
                # se è minore di R_C, salvo la distanza e controllo se è 1° o 2° vicino
                else:
                    d_val = np.sqrt(rij_squared)
                    # se è entro R_P, lo segno come 1° vicino
                    if rij_squared <= R_P_squared:
                        is_neigh = True
                    # se è entro R_C, lo segno come 2° vicino
                    else:
                        is_second = True
                        
            # aggiorno le matrici di vicini e distanze
            neighbour[i, j] = neighbour[j, i] = is_neigh
            second[i, j] = second[j, i] = is_second
            if d_val != np.inf:
                new_dist[i, j] = new_dist[j, i] = d_val
                    
    for k in range(N):
        new_dist[k, k] = 0.0
    return neighbour, second, new_dist

class CrystalStructure:
    """
    Classe per rappresentare una struttura cristallina.
    Attributi:
    - vec_x, vec_y, vec_z: coordinate degli atomi
    - N_atoms: numero totale di atomi
    Metodi:
    - __init__: inizializza gli attributi della classe
    - empty: crea un cristallo vuoto con un numero specificato di atomi
    - from_file: legge un file .txt e restituisce le coordinate degli atomi
    - find_neighbours: trova i primi vicini di ogni atomo in base a una distanza di taglio R_C
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
        self.pcb = np.full(3, 166.416) # *10 così di default non considera la periodicità
        
    @classmethod
    def from_file(cls, filename):
        """Crea un Crystal leggendo da file."""
        data = np.loadtxt(filename)
        return cls(data)
    # Uso:
    # crystal = Crystal.from_file('data.txt')
    
    @classmethod
    def empty(cls, n_atoms):
        """Constructor alternativo: crea Crystal vuoto"""
        zeros = np.zeros(n_atoms)
        return cls(zeros, zeros, zeros)
    
    # viste sulle posizioni
    @property
    def vec_x(self):
        return self.positions[:, 0]

    @property
    def vec_y(self):
        return self.positions[:, 1]

    @property
    def vec_z(self):
        return self.positions[:, 2]
    
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
    # METODI 
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
    
    def copy(self):
        return copy.deepcopy(self)
    
    def set_R_C(self, R_C):
        """
        Imposta la distanza di taglio R_C usata per trovare i vicini.
        """
        self.R_C = R_C
        
    def set_R_P(self, R_P):
        """
        Imposta il punto di giunzione polinomiale R_P; 
        divide primi e secondi vicini.
        """
        self.R_P = R_P
    
    def set_R_V(self, R_V):
        """
        Imposta la grandezza della Verlet cage R_V.
        Consigliato: R_V = R_C + 0.5 Å.
        Non può essere minore di R_C, altrimenti solleva un ValueError.
        """
        if R_V < self.R_C:
            raise ValueError(f"R_V ({R_V}) non può essere minore di R_C ({self.R_C}).")
        self.R_V = R_V
        
    def set_pbc(self, pcb):
        """
        Imposta la dimensione della cella per la condizione di periodicità al contorno.
        Deve essere un array-like di 3 elementi.
        Deve essere consistente con R_C.
        
        pcb: array-like di 3 elementi con le dimensioni della cella in Å.
        """
        # controllo lunghezza
        if len(pcb) != 3:
            raise ValueError("pcb deve essere un array-like di 3 elementi.")
        # controllo consistenza con R_C
        if self.R_C >= 0.5 * min(pcb):
            raise ValueError("⚠️ Attenzione: R_C deve essere minore della metà della dimensione della cella.")
        
        self.pcb = np.asarray(pcb, dtype=np.float64)
        
    def add_atom(self, position):
        """
        Aggiunge un atomo alla struttura cristallina.
        position: array-like di 3 elementi con le coordinate dell'atomo.
        """
        position = np.asarray(position, dtype=np.float64)
        if position.shape != (3,):
            raise ValueError("La posizione deve essere un array-like di 3 elementi.")
        self.positions = np.vstack([self.positions, position])
        self.N_atoms += 1
        self.neighbours_computed = False  # i vicini devono essere ricalcolati
    
    def find_neighbours(self):
        '''
        Ricalcola OGNI distanza e trova primi e secondi vicini per ogni atomo in base a R_C e R_P.
        Gli atomi entro una distanza di Verlet (R_C < r < R_V) sono aggiunti alla matrice dei secondi vicini,
        ma la loro distanza non è salvata in distance_matrix.
        '''
        positions = np.asarray(self.positions, dtype=np.float64)
        neighbour, second, dist = _find_neighbour_masks_kernel(positions, self.R_P, self.R_C, self.R_V, self.pcb)
        # aggiorno le matrici dei vicini e delle distanze
        self.neighbour_matrix = neighbour
        self.second_neighbour_matrix = second
        self.distance_matrix = dist
        # eseguito il calcolo, si attiva un flag e si salvano le posizioni di riferimento
        self.neighbours_computed = True
        self.reference_positions = np.copy(self.positions)
        
    def only_neighbours_distance(self):
        '''
        Ricalcola solamente le distanze dai vicini già noti, senza aggiornarli.
        Se le posizioni sono cambiate troppo, neighbours_computed è posto False ed
        è necessario rieseguire find_neighbours().
        '''
        positions = np.asarray(self.positions, dtype=np.float64)
        neighbour, second, dist = _only_neighbours_distance_kernel(positions, self.R_P, self.R_C, self.R_V, self.neighbour_matrix, self.distance_matrix, self.pcb)
        # aggiorno le matrici dei vicini e delle distanze
        self.neighbour_matrix = neighbour
        self.second_neighbour_matrix = second
        self.distance_matrix = dist
        
    def print_neighbours(self, index=None):
        """
        Stampa gli indici dei primi vicini per ogni atomo o di uno nello specifico.
        """
        if not self.neighbours_computed:
            print("Prima di chiamare questo metodo, esegui find_neighbours(R_C)")
            return None
        
        if index is not None:
            # basta dire quali sono i True sulla riga index di neighbour_matrix:
            print(f"Vicini dell'atomo {index}: {np.where(self.neighbour_matrix[index])[0]}")
            return None
        else:
            print("Indici dei vicini per ogni atomo:")
            for i in range(self.N_atoms):
                print(f"Atomo {i}, n_neigh={np.sum(self.neighbour_matrix[i])}: {np.where(self.neighbour_matrix[i])[0]}")
                
    def print_second_neighbours(self, index=None):
        if not self.neighbours_computed:
            print("Prima di chiamare questo metodo, esegui find_neighbours(R_C)")
            return None 
        if np.isfinite(self.R_P) is False:
            print("R_P non definito; assegnarlo con CrystalStructure.set_R_P ed eseguire find_neighbours().")
            return None
        if index is not None:
            print(f"Secondi vicini dell'atomo {index}: {np.where(self.second_neighbour_matrix[index])[0]}")
            return None
        else:
            print("Indici dei secondi vicini per ogni atomo:")
            for i in range(self.N_atoms):
                print(f"Atomo {i}, n_2nd_neigh={np.sum(self.second_neighbour_matrix[i])}: {np.where(self.second_neighbour_matrix[i])[0]}")