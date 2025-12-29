import numpy as np
from numba import jit
#Lx = Ly = Lz =  16.6416

# restituisce le posizioni degli atomi
def read_pos(filename):
        pos = np.loadtxt(filename)
        natoms= pos.shape[0] 
        return pos, natoms
    
def pbc(pos, Lx, Ly, Lz):
    if Lx > 0:
        pos[:, 0] -= Lx*np.floor(pos[:,0]/Lx)
    if Ly > 0:
        pos[:, 1] -= Ly*np.floor(pos[:,1]/Ly)
    return pos

#if Lz > 0:
    #pos[:, 2] -= Lz*np.floor(pos[:,2]/Lz)        
#----------------------------------------------------------------
#               PRIMI VICINI E DISTANZE
#----------------------------------------------------------------
#restituisce which nbrs nbrs e distanze coi primi vicini

@jit(nopython = True)
def nbrs_and_distances(pos,natoms, rc):
    nbrs = np.zeros(natoms, dtype=np.int64)
    whichnbr = np.zeros((natoms, natoms), dtype=np.int64)
    d2vec = np.zeros((natoms, natoms), dtype = np.float64)
    for i in range(natoms):
        for j in range(natoms):
            if i==j:
                continue
            
            dx = pos[i,0]-pos[j,0]
            dy = pos[i,1]-pos[j,1]
            dz = pos[i,2]-pos[j,2]
            d2 = (dx*dx) + (dy*dy) + (dz*dz)
            d2 = d2**0.5
            if d2 <= rc:
                nbrs[i] +=1
                whichnbr[i, nbrs[i]-1] = j
                d2vec[i, j] = d2 #distanze coi primi vicini
    return nbrs, whichnbr, d2vec

#caso in cui si sta considerando la singola cella, pongo condizioni periodiche al contorno
@jit(cache=False)
def nbrs_and_distances_cell(pos, natoms, rc, Lx, Ly, Lz):
    nbrs = np.zeros(natoms, dtype=np.int64)
    whichnbr = np.zeros((natoms, natoms), dtype=np.int64)
    d2vec = np.zeros((natoms, natoms), dtype = np.float64)
    for i in range(natoms):
        for j in range(natoms):
            if i == j:
                continue
            
            #x------------------
            dx = pos[i,0]-pos[j,0]
            dx = dx - Lx*np.round(dx/Lx)
                
            #y------------------    
            dy = pos[i,1]-pos[j,1]
            dy = dy - Ly*np.round(dy/Ly)
            
            #z--------------    
            dz = pos[i,2]-pos[j,2]
            dz = dz - Lz*np.round(dz/Lz)

            #d--------------------
            d2 = (dx*dx) + (dy*dy) + (dz*dz)
            d2 = d2**0.5
            
            if d2 <= rc:
                nbrs[i] +=1
                whichnbr[i, nbrs[i]-1] = j
                d2vec[i, j] = d2 #distanze coi primi vicini
                    
    return nbrs, whichnbr, d2vec
#----------------------------------------------------------------
#               ENERGIA POTENZIALE
#----------------------------------------------------------------
#restituisce phi per i e j
@jit(cache=False)
def phi_ij(d2: float) -> float:
    eps =  0.345
    sig = 2.644
    phi_ij = 4*eps*((sig/d2)**12 - (sig/d2)**6)
    return phi_ij

#restituisce l'energia potenziale totale

@jit(cache=False)
def E_pot_tot(d2vec, natoms, whichnbr, nbrs, rp, coeffs):
    def poly(coeffs, d2):    
        P7 = coeffs[0] + coeffs[1]*d2 + coeffs[2]*(d2**2) + coeffs[3]*(d2**3) + coeffs[4]*(d2**4) + coeffs[5]*(d2**5) + coeffs[6]*(d2**6) + coeffs[7]*(d2**7)
        return P7
    E_pot = 0
    for i in range(natoms):
        for j in range(nbrs[i]):
            k = whichnbr[i,j]
            if d2vec[i,k] <= rp:
                E_pot += phi_ij(d2vec[i,k])
            else:
                E_pot += poly(coeffs, d2vec[i,k])
    E_pot = E_pot/2
    return E_pot


#----------------------------------------------------------------
#               FORZE
#----------------------------------------------------------------
# restituisce le tre componenti delle forze

@jit(cache=False)
def forces(natoms, nbrs, whichnbr, pos, d2vec, rp, coeffs, Lx, Ly, Lz):
    eps =  0.345
    sig = 2.644
    def poly_der(coeffs, d2):
        P7_prime = coeffs[1] + 2*coeffs[2]*(d2) + 3*coeffs[3]*(d2**2) + 4*coeffs[4]*(d2**3) + 5*coeffs[5]*(d2**4) + 6*coeffs[6]*(d2**5) + 7*coeffs[7]*(d2**6)
        return P7_prime
    
    force = np.zeros((natoms, 3), dtype = np.float64)
    
    for i in range(natoms):
        for j in range(nbrs[i]):
            k = whichnbr[i,j]
            
            dx = pos[i,0]-pos[k,0]
            dy = pos[i,1]-pos[k,1]
            dz = pos[i,2]-pos[k,2]
            
            dx = dx - Lx*np.round(dx/Lx)
            dy = dy - Ly*np.round(dy/Ly)
            dz = dz - Lz*np.round(dz/Lz)

            if d2vec[i,k] <= rp:
                force[i, 0]+= 24*eps*(sig**6)*(d2vec[i,k]**(-8))*(dx)*((2*sig**6)*(d2vec[i,k]**(-6))-1)
                force[i, 1]+= 24*eps*(sig**6)*(d2vec[i,k]**(-8))*(dy)*((2*sig**6)*(d2vec[i,k]**(-6))-1)
                force[i, 2]+= 24*eps*(sig**6)*(d2vec[i,k]**(-8))*(dz)*((2*sig**6)*(d2vec[i,k]**(-6))-1)
            else:
                force[i,0]+= -poly_der(coeffs, d2vec[i,k])*((dx)/d2vec[i,k])
                force[i,1]+= -poly_der(coeffs, d2vec[i,k])*((dy)/d2vec[i,k])
                force[i,2]+= -poly_der(coeffs, d2vec[i,k])*((dz)/d2vec[i,k])
    return force

def forces_old(natoms, nbrs, whichnbr, pos, d2vec, dx, dy, dz, rp, pol, eps =  0.345, sig = 2.644):
    force = np.zeros((natoms, 3))
    for i in range(natoms):
        for j in range(nbrs[i]):
            k = whichnbr[i,j]
            if d2vec[i,k] <= rp:
                force[i, 0]+= 24*eps*(sig**6)*(d2vec[i,k]**(-8))*dx[i,k]*((2*sig**6)*(d2vec[i,k]**(-6))-1)
                force[i, 1]+= 24*eps*(sig**6)*(d2vec[i,k]**(-8))*dy[i,k]*((2*sig**6)*(d2vec[i,k]**(-6))-1)
                force[i, 2]+= 24*eps*(sig**6)*(d2vec[i,k]**(-8))*dz[i,k]*((2*sig**6)*(d2vec[i,k]**(-6))-1)
            else:
                force[i,0]+= -pol.poly_der(d2vec[i,k])*(dx[i,k]/d2vec[i,k])
                force[i,1]+= -pol.poly_der(d2vec[i,k])*(dy[i,k]/d2vec[i,k])
                force[i,2]+= -pol.poly_der(d2vec[i,k])*(dz[i,k]/d2vec[i,k])
    return force

#----------------------------------------------------------------
#               VELOCITA'
#----------------------------------------------------------------
def random_velocities(C, natoms, myseed, T, m):
    k_b = 1/11603
    np.random.seed(myseed)
    vel = np.zeros((natoms, 3))
    vel_x_tot = 0 
    vel_y_tot = 0 
    vel_z_tot = 0
    for i in range(natoms):
        vel[i,0] = np.random.uniform(-C, C)
        vel[i,1] = np.random.uniform(-C, C)
        vel[i,2] = np.random.uniform(-C, C)
    # velocità media
    for i in range(natoms):
        vel_x_tot += vel[i, 0]
        vel_y_tot += vel[i, 1]
        vel_z_tot += vel[i, 2]
        
    vel_x_media = vel_x_tot/natoms
    vel_y_media = vel_y_tot/natoms
    vel_z_media = vel_z_tot/natoms
    
    #correzione alle velocità
    for i in range(natoms):
        vel[i,0] = vel[i,0] - vel_x_media
        vel[i,1] = vel[i,1] - vel_y_media
        vel[i,2] = vel[i,2] - vel_z_media
        
    E_k_tot = kin_energy(vel, m, natoms)
    T_prime = (2*E_k_tot)/(3*natoms*k_b)
    
    for i in range(natoms):
        vel[i,0] = vel[i,0]*(np.sqrt(T/T_prime))
        vel[i,1] = vel[i,1]*(np.sqrt(T/T_prime))
        vel[i,2] = vel[i,2]*(np.sqrt(T/T_prime))
        
    return vel

#----------------------------------------------------------------
#               ENERGIA CINETICA
#----------------------------------------------------------------
def ek_single_atom(v_atomo, mass):
    E_K = 0.5*mass*(v_atomo**2)
    return E_K
    
def kin_energy(vel, mass, natoms):
    E_k_tot = 0
    for i in range(natoms):
        v_atomo = (vel[i,0]**2+vel[i,1]**2+vel[i,2]**2)**0.5
        E_k_tot += ek_single_atom(v_atomo, mass)    
    return E_k_tot

#----------------------------------------------------------------
#               STEEPEST-DESCENT ALGORITHM
#----------------------------------------------------------------
@jit(nopython = True)
def stepest_descent_cell(pos, natoms, coeffs, rc, rp, Csteep, ftoll, maxstep, Lx, Ly, Lz):
    max_force_vec = np.zeros(maxstep, dtype = np.float64)
    e_pot_tot_vec = np.zeros(maxstep, dtype = np.float64) 
    step_convergente = 0
    for i in range(maxstep):
        print(pos)
        nbrs, whichnbr, d2vec = nbrs_and_distances_cell(pos, natoms, rc, Lx, Ly, Lz)
        forze = forces(natoms, nbrs, whichnbr, pos, d2vec, rp, coeffs)
        
        max_force = np.max(np.sqrt(forze[:,0]**2+forze[:,1]**2+forze[:,2]**2))
        max_force_vec[i] =max_force
        epot = E_pot_tot(d2vec, natoms, whichnbr, nbrs, rp, coeffs)
        e_pot_tot_vec[i] = epot 
        
        
        if max_force < ftoll:
            print(f"abbiamo raggiunto la convergenza in {i} steps")
            step_convergente = i
            break # ha raggiunto il minimo
        else:
            pos = pos + (Csteep*forze)

    return pos, max_force_vec, e_pot_tot_vec, step_convergente
    
@jit(nopython = True)
def stepest_descent(pos, natoms, coeffs, rc, rp, Csteep, ftoll, maxstep):
    max_force_vec = np.zeros(maxstep, dtype = np.float64)
    e_pot_tot_vec = np.zeros(maxstep, dtype = np.float64) 
    step_convergente = 0
    for i in range(maxstep):
        nbrs, whichnbr, d2vec = nbrs_and_distances(pos, natoms, rc)
        forze = forces(natoms, nbrs, whichnbr, pos, d2vec, rp, coeffs)
        
        max_force = np.max(np.sqrt(forze[:,0]**2+forze[:,1]**2+forze[:,2]**2))
        max_force_vec[i] =max_force
        epot = E_pot_tot(d2vec, natoms, whichnbr, nbrs, rp, coeffs)
        e_pot_tot_vec[i] = epot 
        
        
        if max_force < ftoll:
            print(f"abbiamo raggiunto la convergenza in {i} steps")
            step_convergente = i
            break # ha raggiunto il minimo
        else:
            pos = pos + (Csteep*forze)

    return pos, max_force_vec, e_pot_tot_vec, step_convergente    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    