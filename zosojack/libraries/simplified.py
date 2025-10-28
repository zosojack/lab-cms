import numpy as np
from libraries.CrystalStructure import CrystalStructure

# Costanti fisiche
k_B = 1/11603 # eV/K

atomic_mass = 108 * 1.66053906660E-27 / 16

def compute_potential(cristallo, sigma=2.644, epsilon=0.345):
            """
            Calcola il potenziale di Lennard-Jones per ogni atomo.
            Il potenziale dipende dai parametri A, B e, indirettamente, dalla distanza di
            taglio R_C, usata per trovare i vicini.
            Restituisce il potenziale totale del sistema.
            """
            def _lennard_jones_ij(r_ij):
                twelve = (sigma/r_ij)**12
                six = (sigma/r_ij)**6
                return 4*epsilon*(twelve - six)
            
            potenziale = 0
            for i, atom_neighbours in enumerate(cristallo.which_neighbour):
                for j, neighbour_index in enumerate(atom_neighbours):
                    if i != neighbour_index:
                        d_ij = cristallo.how_distant[i][j]
                        potenziale += _lennard_jones_ij(d_ij)
                
            cristallo.potential = potenziale / 2
            return cristallo.potential

def compute_forces_matrix(cristallo, sigma=2.644, epsilon=0.345):
        """
        Versione che utilizza matrici numpy.
        Calcola le forze sugli atomi dovute al potenziale di Lennard-Jones (F = -∇V).
        Restituisce una matrice Nx3 di forze per ogni atomo.
        """
        
        def addendo_derivata_lennard_jones(q_i, q_k, r_ik):
            return 1/(r_ik**8) * ( (2*sigma**6)/(r_ik**6) - 1 ) * (q_i - q_k)

        mat_forza = np.zeros((cristallo.N_atoms, 3)) # ciascuna riga è un atomo, ognuna con tre colonne
        for i, atom_neighbours in enumerate(cristallo.which_neighbour):
            for j, neighbour_index in enumerate(atom_neighbours):
                d_ij = cristallo.how_distant[i][j]
                mat_forza[i, 0] += addendo_derivata_lennard_jones(cristallo.vec_x[i], cristallo.vec_x[neighbour_index], d_ij)
                mat_forza[i, 1] += addendo_derivata_lennard_jones(cristallo.vec_y[i], cristallo.vec_y[neighbour_index], d_ij)
                mat_forza[i, 2] += addendo_derivata_lennard_jones(cristallo.vec_z[i], cristallo.vec_z[neighbour_index], d_ij)
        mat_forza *= 24 * epsilon * sigma**6

        return mat_forza


def random_velocities(cristallo, temp_ini):
    """
    Inizializza le velocità casuali secondo una distribuzione uniforme
    tra -C e C, con C = √(3*k_B*T/m). È praticamente K(0). 
    E_tot = K(0) + V(0) è stabilita dalla scelta di T.
    Restituisce una matrice Nx3 di velocità.
    """
    atomic_mass = 108 * 1.66053906660E-27 / 16
    
    N_atoms = cristallo.N_atoms 
    C = np.sqrt((3*k_B*temp_ini) / atomic_mass)
    vel = np.zeros((N_atoms, 3))
    
    for i in range(N_atoms):
        vel[i,0] = np.random.uniform(-C, C)
        vel[i,1] = np.random.uniform(-C, C)
        vel[i,2] = np.random.uniform(-C, C)
        
    vel_x_media = np.mean(vel[:,0])
    vel_y_media = np.mean(vel[:,1])     
    vel_z_media = np.mean(vel[:,2])     
    
    E_K = 0
    # 2 CORREZIONI
    # - prima correzione del drift sottraendo la media
    # - poi rescale dovuto alla modifica di T con fattore correttivo √(T_ini / T_prime)
    for i in range(N_atoms):
        # prima correzione del drift
        vel[i,0] = vel[i,0] - vel_x_media
        vel[i,1] = vel[i,1] - vel_y_media
        vel[i,2] = vel[i,2] - vel_z_media
        # calcolo energia cinetica
        E_K += 0.5*atomic_mass*(vel[i,0]**2 + vel[i,1]**2 + vel[i,2]**2)
    T_prime = (2/3) * (E_K / (N_atoms * k_B)) # temperatura effettiva dopo la prima correzione
    
    velocities = vel * np.sqrt(temp_ini / T_prime) 
    kinetic_E = E_K
    
    return velocities, kinetic_E
    
def update_positions(cristallo, velocities, forze, dt):
    """
    Aggiorna le posizioni degli atomi usando il metodo di Verlet.
    x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt^2
    """
    acc = forze / atomic_mass  # accelerazioni
    cristallo.positions[:, 0] += velocities[:, 0] * dt + 0.5 * acc[:, 0] * dt**2
    cristallo.positions[:, 1] += velocities[:, 1] * dt + 0.5 * acc[:, 1] * dt**2
    cristallo.positions[:, 2] += velocities[:, 2] * dt + 0.5 * acc[:, 2] * dt**2

def update_velocities(cristallo, velocità, forza_prima, forza_dopo, dt):
    """
    Aggiorna le velocità degli atomi usando il metodo di Verlet.
    Aggiorna anche il valore dell'energia cinetica.
    v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt
    """
    
    acc_old = forza_prima / atomic_mass
    acc_new = forza_dopo / atomic_mass
    E_K = 0
    
    for i in range(cristallo.N_atoms):
        # Aggiorno velocità
        velocità[i, 0] += 0.5 * (acc_old[i, 0] + acc_new[i, 0]) * dt
        velocità[i, 1] += 0.5 * (acc_old[i, 1] + acc_new[i, 1]) * dt
        velocità[i, 2] += 0.5 * (acc_old[i, 2] + acc_new[i, 2]) * dt
        # Aggiorno energia cinetica
        E_K += 0.5*atomic_mass*(velocità[i, 0]**2 + velocità[i, 1]**2 + velocità[i, 2]**2)

    return velocità, E_K


