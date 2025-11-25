# CrystalStructure.py 
import numpy as np
from numba import njit

@njit(cache=True)
def _neighbour_masks_kernel(pos, R_P, R_C):
    N = pos.shape[0]
    neighbour = np.zeros((N, N), dtype=np.bool_)
    second = np.zeros((N, N), dtype=np.bool_)
    dist = np.full((N, N), np.inf)
    for i in range(N):
        xi, yi, zi = pos[i, 0], pos[i, 1], pos[i, 2]
        for j in range(i + 1, N):
            dx = xi - pos[j, 0]
            dy = yi - pos[j, 1]
            dz = zi - pos[j, 2]
            rij = np.sqrt(dx*dx + dy*dy + dz*dz)
            if rij <= R_C:
                dist[i, j] = dist[j, i] = rij
                if rij <= R_P:
                    neighbour[i, j] = neighbour[j, i] = True
                else:
                    second[i, j] = second[j, i] = True
    for k in range(N):
        dist[k, k] = 0.0
    return neighbour, second, dist

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
        self.N_neighbours = None  # numero di primi vicini per ogni atomo
        self.which_neighbour = None  # indici dei primi vicini per ogni atomo
        self.how_distant = None  # distanze tra i primi vicini per ogni atomo
        self.R_C = np.inf  # distanza di taglio usata per trovare i vicini
        self.R_P = np.inf  # punto di giunzione polinomiale
        
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
    
    def find_neighbours_numba(self):
        pos = np.asarray(self.positions, dtype=np.float64)
        neighbour, second, dist = _neighbour_masks_kernel(pos, self.R_P, self.R_C)
        self.neighbour_matrix = neighbour
        self.second_neighbour_matrix = second
        self.distance_matrix = dist
    

    def print_neighbours(self, index=None):
        """
        Stampa l'attributo which_neighbour che contiene gli indici dei primi vicini per ogni atomo.
        """
        if self.which_neighbour is None:
            print("Prima di chiamare questo metodo, esegui find_neighbours(R_C)")
            return None
        
        if index is not None:
            print(f"Vicini dell'atomo {index}: {self.which_neighbour[index]}")
            return None
        else:
            print("Indici dei vicini per ogni atomo:")
            for i, neighbours in enumerate(self.which_neighbour):
                print(f"Atomo {i}, n_neigh={self.N_neighbours[i]}: {neighbours}")