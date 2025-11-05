# CrystalDynamics.py

import numpy as np
import os
from pathlib import Path
from libraries.CrystalStructure import CrystalStructure
from libraries.CrystalPotential import CrystalPotential

# Costanti fisiche
k_B = 1/11603 # eV/K
N_A = 1.66053906660E-27 # numero di Avogadro

'''
Deve avere:
- funzione che setta le velocità iniziali secondo dist uniforme ok
- funzione che calcola l'energia cinetica totale 
- funzione che aggiorna le posizioni e velocità secondo Verlet ok

Prende in input un oggetto CrystalStructure, che ha la struttura iniziale
devo pensare poi se ha senso aggiornarne o no gli attributi. 
-> sì dai. facciamo in modo che CrystalDynamics modifichi l'oggetto CrystalStructure
passato in input. Agirà sulla matrice delle coordinate.


Ci sta spostare la funzione compute_forces qui, così come anche potential_energy.
Perché andranno chiamate più volte durante i loop di dinamica.
'''

class CrystalDynamics:
    
    def __init__(self, 
                 crystal: CrystalStructure, 
                 atomic_mass: float = 108, 
                 dt: float = 1E-15,
                 temp_ini: float = 20):
        # Oggetti 
        self.crystal  = crystal  # oggetto CrystalStructure
        # Parametri della simulazione - FIXME: questi meglio come argomento di run_dynamics?
        self.atomic_mass = atomic_mass * 1.66E-27 / 16  # massa degli atomi in eV s^2/Angstrom^2
        self.dt = dt      # passo temporale
        self.temp_ini = temp_ini  # temperatura iniziale della simulazione
        # Variabili di stato
        self.old_force = None  # forze sugli atomi tempo t
        self.new_force = None  # forze sugli atomi tempo t+dt
        self.velocities = None  # velocità
        self.temp_final = None  # temperatura finale della simulazione
        self.kinetic_E = None
        self.potential_E = None
        
          
    '''
    DINAMICA
    0) Inizio
    - Leggo posizioni iniziali
    - Stabilisco vicini
    - Calcolo potenziale iniziale
    - Calcolo forze iniziali
    - Inizializzo velocità con random
    1) Dinamica
    - Aggiorno posizioni
    - Calcolo forze nuove
    - Aggiorno velocità
    - Calcolo energie

    '''

    def set_seed(self, myseed):
        """
        Imposta il seed per la generazione di numeri casuali.
        """
        np.random.seed(myseed)

    def _random_velocities(self):
        """
        Inizializza le velocità casuali secondo una distribuzione uniforme
        tra -C e C, con C = √(3*k_B*T/m). È praticamente K(0). 
        E_tot = K(0) + V(0) è stabilita dalla scelta di T.
        Restituisce una matrice Nx3 di velocità.
        """
        N_atoms = self.crystal.N_atoms 
        C = np.sqrt((3*k_B*self.temp_ini) / self.atomic_mass)
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
            E_K += 0.5*self.atomic_mass*(vel[i,0]**2 + vel[i,1]**2 + vel[i,2]**2)
            
        T_prime = (2/3) * (E_K / (N_atoms * k_B)) # temperatura effettiva dopo la prima correzione
        
        self.velocities = vel * np.sqrt(self.temp_ini / T_prime) 
        self.kinetic_E = E_K
        
    def _update_positions(self):
        """
        Aggiorna le posizioni degli atomi usando il metodo di Verlet.
        x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt^2
        """
        acc = self.old_force / self.atomic_mass  # accelerazioni
        self.crystal.positions += self.velocities * self.dt + 0.5 * acc * (self.dt**2)
        
    def _update_neighbours(self):
        """
        Aggiorna la lista dei vicini dopo aver modificato le posizioni.
        """
        self.crystal.find_neighbours()
        
    def _update_forces(self):
        """
        Calcola le nuove forze sugli atomi.
        """
        if self.new_force is not None:
            self.old_force = self.new_force
        self.new_force = CrystalPotential(self.crystal).compute_forces_matrix()
    
    def _update_velocities(self):
        """
        Aggiorna le velocità degli atomi usando il metodo di Verlet.
        Aggiorna anche il valore dell'energia cinetica.
        v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt
        """
        acc_old = self.old_force / self.atomic_mass
        acc_new = self.new_force / self.atomic_mass
        E_K = 0
        vel = self.velocities.copy()
        
        for i in range(self.crystal.N_atoms):
            # Aggiorno velocità
            vel[i, 0] += 0.5 * (acc_old[i, 0] + acc_new[i, 0]) * self.dt
            vel[i, 1] += 0.5 * (acc_old[i, 1] + acc_new[i, 1]) * self.dt
            vel[i, 2] += 0.5 * (acc_old[i, 2] + acc_new[i, 2]) * self.dt
            # Aggiorno energia cinetica
            E_K += 0.5*self.atomic_mass*(vel[i, 0]**2 + vel[i, 1]**2 + vel[i, 2]**2)
            
        self.velocities = vel  # aggiorna la velocità corrente
        self.kinetic_E = E_K # aggiorna energia cinetica corrente

    def _temperature(self):
        '''
        Calcola la temperatura del sistema al tempo t.
        T(t) = 2/3 * K(t) / (N_atoms * k_B)
        '''
        if self.kinetic_E is not None:
            return (2/3) * (self.kinetic_E / (self.crystal.N_atoms * k_B))
        else:
            return self.temp_ini


    def _output_state(self, filename, step, E_tot, E_pot, E_kin, temp):
        """
        Salva lo stato attuale della simulazione in un file. 
        L'i-esima riga corrisponde allo stato istantaneo all'i-esimo step.
        Formato: 
        tempo, energia totale, energia potenziale, energia cinetica, temperatura
        """
        with open(filename, 'w') as f:
            f.write(f"{step*self.dt} \t {E_tot} \t {E_pot} \t {E_kin} \t {temp}")
        
    
    def _output_positions(self, foldername, step, n_steps):
        """
        Salva le posizioni attuali di tutti gli atomi della simulazione in un file.
        Ciascun file corrisponde a uno step.
        L'i-esima riga corrisponde alle posizioni istantanee dell'i-esimo-2 atomo.
        """
        width = len(str(n_steps))          
        step_str = f"{step:0{width}d}"  
        
        filename = foldername \
            / (f"t={self.dt}~" + step_str + f"_{n_steps}.xyz")
        header = f"{self.crystal.N_atoms}\n time={step*self.dt}"
        np.savetxt(filename, self.crystal.positions, header=header, comments='')
        
        
    

    # FIXME: pensare se mettere temperatura e dt come argomenti di run_dynamics
    def run_dynamics(self, n_steps, output=False, debug=False):
        """
        Esegue la simulazione di dinamica molecolare per n_steps.
        """    
        # Predispone la cartella di output, se necessaria
        # TODO: potrebbe essere necessario aggiungere il seed 
        if output:
            out_dir = Path(f"output/dynamics/steps{n_steps}~dt{self.dt}~T{self.temp_ini}~Ag~{self.crystal.N_atoms}")
            out_dir.mkdir(parents=True, exist_ok=True) # crea se manca
            state_file = out_dir / "energy.txt"
            
            n_Print = 10
        
        # Controlla che i vicini siano stati calcolati
        if self.crystal.which_neighbour is None:
            print(f"⚠️  Vicini non calcolati in precedenza. Calcolo con R_C={self.crystal.R_C}.")
            self.crystal.find_neighbours()
        
        # Inizializza velocità casuali
        if self.velocities is None:
            self._random_velocities()
        
        # Calcola forze iniziali
        if self.old_force is None:
            self.old_force = CrystalPotential(self.crystal).compute_forces_matrix()
            
        # Oggetti che raccolgono i metadati (TODO: forse meglio np.array)
        meta_E_tot = []
        meta_E_K = []
        meta_T = []
        
        # LOOP PRINCIPALE
        # qui andrà scartata la fase iniziale. 
        # finché la tempreratura non si stabilizza a T(t) ~ T(0)/2  
        
        for step in range(n_steps):
            # Aggiorna posizioni
            self._update_positions()
            # Aggiorna vicini
            self._update_neighbours()
            # Calcola nuove forze
            self._update_forces()
            # Aggiorna velocità
            self._update_velocities()
            # calcola energie e temperatura parziale
            potential_energy = CrystalPotential(self.crystal).compute_potential()
            temp = self._temperature()
            
            E_tot_now = potential_energy + self.kinetic_E
            meta_E_tot.append(E_tot_now)
            meta_E_K.append(self.kinetic_E)
            meta_T.append(temp)
            
            if debug:
                print(f"step {step+1}/{n_steps}: E_tot={E_tot_now:.3f} eV, V={potential_energy:.3f} eV, K={self.kinetic_E:.3f} eV, T={temp:.1f} K")
                
            
            if output:
                if step % n_Print == 0:
                    self._output_state(state_file,
                                       step,
                                       E_tot_now,
                                       potential_energy,
                                       meta_E_K,
                                       meta_T)
                    self._output_positions(out_dir, step, n_steps)

        return meta_E_tot, meta_E_K, meta_T
        