import numpy as np

class Crystal:
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
    def __init__(self, vec_x, vec_y, vec_z):
        self.vec_x = vec_x  # coordinate x degli atomi
        self.vec_y = vec_y  # coordinate y degli atomi
        self.vec_z = vec_z  # coordinate z degli atomi
        self.N_atoms = len(vec_x)  # numero totale di atomi
        self.N_neighbours = None  # numero di primi vicini per ogni atomo
        self.which_neighbour = None  # indici dei primi vicini per ogni atomo
        self.how_distant = None  # distanze tra i primi vicini per ogni atomo
        self.potential = None  # potenziale calcolato 
        
    @classmethod
    def from_file(cls, filename):
        """Crea un Crystal leggendo da file."""
        data = np.loadtxt(filename)
        vec_x = data[:, 0] 
        vec_y = data[:, 1]  
        vec_z = data[:, 2]
        return cls(vec_x, vec_y, vec_z)
    # Uso:
    # crystal = Crystal.from_file('data.txt')
    
    @classmethod
    def empty(cls, n_atoms):
        """Constructor alternativo: crea Crystal vuoto"""
        zeros = np.zeros(n_atoms)
        return cls(zeros, zeros, zeros)

    def find_neighbours(self, R_C):
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
                    
                    if d_ij <= R_C:
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
        
    def calculate_potential(self):
        """
        Calcola il potenziale di Lennard-Jones per ogni atomo.
        Il potenziale dipende dai parametri A, B e dalla distanza di taglio R_C.
        """
        def lennard_jones_ij(r_ij, sigma=2.644, epsilon=0.345):
            twelve = (sigma/r_ij)**12
            six = (sigma/r_ij)**6
            return 4*epsilon*(twelve - six)
        
        potenziale = 0
        for i, atom_neighbours in enumerate(self.which_neighbour):
            for j, neighbour_index in enumerate(atom_neighbours):
                if i != neighbour_index:
                    d_ij = self.how_distant[i][j]
                    potenziale += lennard_jones_ij(d_ij)
            
        self.potential = potenziale / 2
        return self.potential
    
    def calculate_forces(self):
        
        sigma = 2.644  # parametro sigma per Lennard-Jones
        epsilon = 0.345  # parametro epsilon per Lennard-Jones
        
        def addendo_derivata_lennard_jones(q_i, q_k, r_ik):
            return 1/(r_ik**8) * ( (2*sigma**6)/(r_ik**6) - 1 ) * (q_i - q_k)
            
        vec_forza = [] # ciascun entrata è un atomo, ciascun atomo è una lista di tre entrate
        for i, atom_neighbours in enumerate(self.which_neighbour):
            forza_x, forza_y, forza_z = 0, 0, 0
            for j, neighbour_index in enumerate(atom_neighbours):
                d_ij = self.how_distant[i][j]
                forza_x += addendo_derivata_lennard_jones(self.vec_x[i], self.vec_x[neighbour_index], d_ij)
                forza_y += addendo_derivata_lennard_jones(self.vec_y[i], self.vec_y[neighbour_index], d_ij)
                forza_z += addendo_derivata_lennard_jones(self.vec_z[i], self.vec_z[neighbour_index], d_ij)
            forza_x *= 24 * epsilon * sigma**6
            forza_y *= 24 * epsilon * sigma**6
            forza_z *= 24 * epsilon * sigma**6
            vec_forza.append([forza_x, forza_y, forza_z])
            
        return vec_forza