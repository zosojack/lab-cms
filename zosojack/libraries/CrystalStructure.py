# CrystalStructure.py 

import numpy as np

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

    # NOTE: questo metodo usa liste, non matrici numpy.
    # utile per capire come funziona, ma find_neighbours_numpy è da preferire.
    def find_neighbours(self):
        """
        Trova i primi vicini di ogni atomo in base a una distanza di taglio R_C.
        Restituisce il numero di vicini, gli indici dei vicini e le distanze.
        """
        
        N_neighbours = []  # lista del numero di primi vicini per ogni atomo
        which_neighbour = []  # lista degli indici dei primi vicini per ogni atomo
        how_distant = []  # lista delle distanze tra i primi vicini per ogni atomo

        for i, (x_i, y_i, z_i) in enumerate(zip(self.vec_x, self.vec_y, self.vec_z)):
            n_neigh_i = 0  # contatore dei primi vicini dell'i-esimo atomo
            list_neigh_i = [] # lista di indici dei primi vicini dell'i-esimo atomo
            list_dist_ij = [] # lista di distanze tra i primi vicini 

            for j, (x_j, y_j, z_j) in enumerate(zip(self.vec_x, self.vec_y, self.vec_z)):
                if i != j:  # non considero la distanza di un atomo da se stesso
                    dx = (x_i - x_j)**2 
                    dy = (y_i - y_j)**2
                    dz = (z_i - z_j)**2
                    d_ij = np.sqrt(dx + dy + dz)

                    if d_ij <= self.R_C:
                        n_neigh_i += 1
                        list_neigh_i.append(j)
                        list_dist_ij.append(d_ij)
                        
            # Ordino i vicini per distanza crescente
            if list_neigh_i:  # se ci sono vicini
                # Creo coppie (distanza, indice) e le ordino per distanza
                paired = list(zip(list_dist_ij, list_neigh_i))
                paired_sorted = sorted(paired, key=lambda x: x[0])
                # Separo nuovamente distanze e indici ordinati
                list_dist_ij = [dist for dist, _ in paired_sorted]
                list_neigh_i = [idx for _, idx in paired_sorted]
            
            # Aggiungo i risultati per l'atomo i-esimo
            N_neighbours.append(n_neigh_i)
            which_neighbour.append(list_neigh_i)
            how_distant.append(list_dist_ij)
            
        self.N_neighbours = N_neighbours
        self.which_neighbour = which_neighbour
        self.how_distant = how_distant
        

    def find_neighbours_matrix(self):
        """
        Versione con matrici numpy di find_neighbours
        - neighbour_matrix: matrice simmetrica di bool che indica se due atomi sono vicini
        - distance_matrix: matrice simmetrica delle distanze tra atomi (inf se non vicini)
        """
        
        n = self.N_atoms
        # Matrice simmetrica: neighbour_matrix[i,j] = True se j è vicino di i
        neighbour_matrix = np.zeros((n, n), dtype=bool) 
        distance_matrix = np.full((n, n), np.inf)  # distanze
        
        for i in range(n):
            for j in range(i+1, n):  # solo triangolo superiore
                dx = self.vec_x[i] - self.vec_x[j]
                dy = self.vec_y[i] - self.vec_y[j] 
                dz = self.vec_z[i] - self.vec_z[j]
                d_ij = np.sqrt(dx**2 + dy**2 + dz**2)

                if d_ij <= self.R_C:
                    neighbour_matrix[i, j] = True
                    neighbour_matrix[j, i] = True  # simmetria
                    distance_matrix[i, j] = d_ij
                    distance_matrix[j, i] = d_ij
        
        self.neighbour_matrix = neighbour_matrix
        self.distance_matrix = distance_matrix 
        

    def find_neighbours_numpy(self):
        """
        Versione con matrici numpy di find_neighbours, ancora più efficiente.
        - neighbour_matrix: matrice simmetrica di bool che indica se due atomi sono vicini
        - distance_matrix: matrice simmetrica delle distanze tra atomi (inf se non vicini)
        """
        N = self.N_atoms
        pos = np.asarray(self.positions)              # (N,3)
        # matrice delle distanze (N,N)
        diffs = pos[:, None, :] - pos[None, :, :]     # (N,N,3)
        dist = np.linalg.norm(diffs, axis=2)
        # mask vicini (escludo la diagonale)
        mask = (dist <= self.R_C) & (~np.eye(N, dtype=bool))
        # risultati
        self.neighbour_matrix = mask
        self.distance_matrix = np.where(mask, dist, np.inf)
        np.fill_diagonal(self.distance_matrix, 0.0)


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