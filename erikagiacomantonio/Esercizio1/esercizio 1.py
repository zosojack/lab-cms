#esercizio 1

# devo capire in che cartella sono perché non mi importa le librerie


# Verlet algorithm
import cristallo as cr 
import numpy as np
from numba import njit
import os
from PolynomialJunction import PolynomialJunction as poly7

Lx = 16.6416
Ly = 16.6416
Lz = 1664
L = np.array([Lx, Ly, Lz])

eps =  0.345
sig = 2.644 
t_ex = 850*2 #la tmeperatura iniziale deve essere tale per cui quella di equilibrio sia 850K
Temp_in = [t_ex] 
Temp_in = np.array(Temp_in)
rc = 4.5
rp = 4.2
coeffs = poly7(epsilon=eps,sigma=sig,R_C=rc,R_P=rp).coeffs_array
k_b = 1/11603
massa_ag = 108*((1.66*(10**(-27)))/16)
dt = [1E-15]
#dt =  np.linspace(1E-15, 10E-15, 8) #1 fs
dt = np.array(dt)
ndt= dt.shape[0]
ntemp= Temp_in.shape[0]
iter_steps = int(300E-12/1E-15)
#int(300E-12/3E-15) #la mia simulazione deve durare 300ps
print(iter_steps) 

myseed = 12316543


filename = "/Users/zosojack/lab-cms/erikagiacomantonio/Esercizio1/fcc100a256_1.txt"
nprint = 100
noutput = 1000

# update positions
@njit(cache=False)
def update_positions(pos: np.ndarray, vel: np.ndarray, forze: np.ndarray, massa_ag: float, dt: float) -> np.ndarray:
    for i in range(pos.shape[0]):
        for j in range(3):
            pos[i,j] = pos[i,j] + vel[i,j]*dt + 0.5*(forze[i,j]/massa_ag)*(dt**2)
    return pos

# update velocities
@njit(cache=False)
def update_velocities(vel: np.ndarray, forces: np.ndarray, forces_new:np.ndarray, massa_ag: float, dt: float) -> np.ndarray:
    for i in range(vel.shape[0]):
        for j in range(3):
            vel[i,j] = vel[i,j] + (forces[i,j] + forces_new[i,j])/(2*massa_ag)*dt
    return vel


#os.makedirs("risultati_esercizio1_2_MD")  #NEL CASO IN CUI SERVISSE CREARE LA CARTELLA TOGLILA DAL COMMENTO
for ntempi in range(ntemp):
    print(f"risultati_esercizio1_4_MD/temperatura iniziale di {Temp_in[ntempi]}K")
    const = ((3*Temp_in[ntempi]*k_b)/massa_ag)**0.5
    pos, natoms = cr.read_pos(filename)
    # applicazione delle periodic boundary conditions
    pos = cr.pbc(pos, Lx, Ly, Lz)
    #nbrs, whichnbr, d2vec = cr.nbrs_and_distances(pos, natoms, rc)
    nbrs, whichnbr, d2vec = cr.nbrs_and_distances_cell(pos, natoms, rc, Lx, Ly, Lz)
    forze = cr.forces(natoms, nbrs, whichnbr, pos, d2vec, rp, coeffs)
    vel = cr.random_velocities(const, natoms, myseed, Temp_in[ntempi], massa_ag)
    epot = cr.E_pot_tot(d2vec, natoms, whichnbr, nbrs, rp, coeffs)
    ektot = cr.kin_energy(vel, massa_ag, natoms)
    etot_vec = np.zeros((iter_steps,2))
    tempo_vec = np.zeros((iter_steps,2))
    print(f"la mia C  è {const:.10e}")
    E_tot = epot + ektot
    print(f"l'energia iniziale totale è {E_tot}")
    print(f"l'energia potenziale totale è {epot}")
    
    os.makedirs(f"erikagiacomantonio/risultati_esercizio1_4_MD/risultati_{Temp_in[ntempi]:.0f}K", exist_ok=True)
    for dti in range(ndt):
        sec = dt[dti]*1E15
        os.makedirs(f"erikagiacomantonio/risultati_esercizio1_4_MD/risultati_{Temp_in[ntempi]:.0f}K/risultati_timestep{sec:.0f}fs", exist_ok=True)
        print(f"----------steps da {dt[dti]:.2e} s")
        #file di testo
        fout = open(f'erikagiacomantonio/risultati_esercizio1_4_MD/risultati_{Temp_in[ntempi]:.0f}K/risultati_timestep{sec:.0f}fs/output_{sec:.0f}fs_{Temp_in[ntempi]:.0f}K.txt', 'w')
        
        fout.write(f'simulazione con {natoms} atomi\n')
        fout.write('time, etot, epot, ekin, temp\n')
        
        print("---------------------------- INIZIO SIMULAZIONE ----------------------------")
        
        for iters in range(iter_steps):
            
            print(f"+++ inizio step {iters} +++") 
            #aggiorno le posizioni
            pos = update_positions(pos, vel, forze, massa_ag, dt[dti])
            
            # ora calcolo le cose nuove con le posizioni nuove
            #nbrs, whichnbr, d2vec = cr.nbrs_and_distances(pos, natoms, rc)
            nbrs, whichnbr, d2vec = cr.nbrs_and_distances_cell(pos, natoms, rc, Lx, Ly, Lz)
            forces_new = cr.forces(natoms, nbrs, whichnbr, pos, d2vec, rp, coeffs)
            
            vel = update_velocities(vel, forze, forces_new, massa_ag, dt[dti])
                    
            ektot = cr.kin_energy(vel, massa_ag, natoms)
            epot = cr.E_pot_tot(d2vec, natoms, whichnbr, nbrs, rp, coeffs)
            E_tot = ektot + epot
            temp = (2*ektot)/(3*natoms*k_b)
            tempo = (dt[dti])*(iters)
            
            if (iters % nprint == 0):
                
                np.savetxt(f'erikagiacomantonio/risultati_esercizio1_4_MD/risultati_{Temp_in[ntempi]:.0f}K/risultati_timestep{sec:.0f}fs/{Temp_in[ntempi]:.0f}K_{sec:.0f}fs_{iters}.xyz', pos, header=f'{natoms}\n time={tempo}', comments = '')
                #print(f"allo step {iters} l'energia totale è {E_tot}")
                
            if (iters % noutput == 0):
                print(f"allo step {iters} l'energia totale è {E_tot}")
                      
            #le forze nuove diventano le forze vecchie
            forze = forces_new
            
            #vettore con un numero di colonne pari al numero di tempi (mi sarebbe servito x fare il graph)
            #etot_vec[iters, dti] = E_tot
            #tempo_vec[iters, dti] = tempo
            
        
            fout.write(f'{tempo} \t {E_tot} \t {epot} \t {ektot} \t {temp}\n')
            
            print(f"+++ fine step {iters} +++") 
            
        print("----------------------------- FINE SIMULAZIONE ----------------------------")
            
    fout.close()

