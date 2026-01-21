# NumbaSubroutines.py
import numpy as np
from numba import njit

# calcoli preliminari: i vicini
@njit
def pbc_corr (xi: float, 
              Li: float
) -> int:
    ''' Correzione condizione periodica al contorno '''
    
    xi -= Li * np.floor(xi/Li)
    return int(xi)

@njit
def neigh_XY (x_i: float,
              y_i: float,
              height: np.ndarray # potrei passare solo una slice della matrice, invece che tutta
) -> int:
    ''' Conta quanti primi vicini possiede l'atomo sulla superficie in posizione X,Y '''
    
    Lx, Ly = height.shape[0], height.shape[1]
    n = 0
    z_c = height[x_i, y_i]
    
    # Vicino 1 (Up)
    if height[x_i, pbc_corr(y_i+1, Ly)] >= z_c: n += 1
    # Vicino 2 (Down)
    if height[x_i, pbc_corr(y_i-1, Ly)] >= z_c: n += 1
    # Vicino 3 (Right)
    if height[pbc_corr(x_i+1, Lx), y_i] >= z_c: n += 1
    # Vicino 4 (Left)
    if height[pbc_corr(x_i-1, Lx), y_i] >= z_c: n += 1
            
    return int(n)

@njit
def count_NN (height: np.ndarray) -> np.ndarray:
    ''' Produce una matrice di interi che contiene il numero di primi vicini per ciascun atomo sulla superficie '''
    
    first_neigh = np.zeros_like(height)
    Lx, Ly = height.shape[0], height.shape[1]
    
    for x in range(Lx):
        for y in range(Ly):
            first_neigh[x,y] = neigh_XY(x, y, height)
    
    return first_neigh

@njit
def _update_NN_deposition (first_neigh: np.ndarray,
                           k_diff: np.ndarray,                # matrice dei rate di diffusione
                           current_k_diff_sum: float,         # somma prima dell'aggiornamento
                           height: np.ndarray,
                           deposition_site: tuple[int, int],
                           params: tuple                      # (J0, J1, T, nu, k_B)
) -> np.ndarray:
    ''' 
    Aggiorna i primi vicini e la matrice dei rate di diffusione dopo un evento di deposizione;
    agisce solamente sui siti coinvolti nell'evento.
    '''
    x0, y0 = deposition_site
    Lx, Ly = height.shape[0], height.shape[1]
    # scorporo i parametri
    J0, J1, T, nu, k_B = params
    
    sites_to_update = [
        (x0, y0), # deposition site
        (x0, pbc_corr(y0+1, Ly)), # Up
        (x0, pbc_corr(y0-1, Ly)), # Down
        (pbc_corr(x0+1, Lx), y0), # Right
        (pbc_corr(x0-1, Lx), y0), # Left
    ]
    
    # è sufficiente aggiornare 5 siti coinvolti (deposition_site + 4 vicini)
    for site in sites_to_update:
        x_dep, y_dep = site
        # aggiorna primo vicino
        first_neigh[x_dep, y_dep] = neigh_XY(x_dep, y_dep, height)
        # aggiorna rate di diffusione
        current_k_diff_sum -= k_diff[x_dep, y_dep] # rimuovo il vecchio contributo
        Eb = -(J0 + first_neigh[x_dep, y_dep] * J1)
        k_diff[x_dep, y_dep] = 4 * nu * np.exp(-Eb / (k_B*T))
        current_k_diff_sum += k_diff[x_dep, y_dep] # aggiungo il nuovo contributo
    
    return first_neigh, k_diff, current_k_diff_sum
    
@njit
def _update_NN_diffusion  (first_neigh: np.ndarray,
                           k_diff: np.ndarray,              # matrice dei rate di diffusione
                           current_k_diff_sum: float,       # somma prima dell'aggiornamento         
                           height: np.ndarray,
                           diffusion_from: tuple[int, int],
                           diffusion_to: tuple[int, int],
                           params: tuple                    # (J0, J1, T, nu, k_B)    
) -> np.ndarray:
    ''' 
    Aggiorna i primi vicini e la matrice dei rate di diffusione dopo un evento di diffusione;
    agisce solamente sui siti coinvolti nell'evento.
    '''
    x_from, y_from = diffusion_from
    x_to, y_to = diffusion_to
    Lx, Ly = height.shape[0], height.shape[1]
    # scorporo i parametri
    J0, J1, T, nu, k_B = params
    
    sites_to_update = [
        (x_from, y_from), # diffusion from site
        (x_from, pbc_corr(y_from+1, Ly)), # Up from
        (x_from, pbc_corr(y_from-1, Ly)), # Down from
        (pbc_corr(x_from+1, Lx), y_from), # Right from
        (pbc_corr(x_from-1, Lx), y_from), # Left from
        (x_to, y_to), # diffusion to site
        (x_to, pbc_corr(y_to+1, Ly)), # Up to
        (x_to, pbc_corr(y_to-1, Ly)), # Down to
        (pbc_corr(x_to+1, Lx), y_to), # Right to
        (pbc_corr(x_to-1, Lx), y_to), # Left to
    ]
    
    # è sufficiente aggiornare i 10 siti coinvolti (diffusion_from + diffusion_to + 8 vicini)
    for site in sites_to_update:
        x_site, y_site = site
        # aggiorna primo vicino
        first_neigh[x_site, y_site] = neigh_XY(x_site, y_site, height)
        # aggiorna rate di diffusione
        current_k_diff_sum -= k_diff[x_site, y_site] # rimuovo il vecchio contributo
        Eb = -(J0 + first_neigh[x_site, y_site] * J1)
        k_diff[x_site, y_site] = 4 * nu * np.exp(-Eb / (k_B*T))
        current_k_diff_sum += k_diff[x_site, y_site] # aggiungo il nuovo contributo
    
    return first_neigh, k_diff, current_k_diff_sum

# selezione evento di diffusione
@njit
def find_move (rho: float,
               k_diff: np.ndarray
) -> tuple:
    ''' 
    Percorre la matrice dei rate di diffusione fino a raggiungere un valore >= rho;
    restituisce le coordinate (x,y) dell'atomo da muovere.
    '''
    Lx, Ly = k_diff.shape[0], k_diff.shape[1]
    par_sum = 0
    
    for x in range(Lx):
        for y in range(Ly):
            par_sum += k_diff[x,y]
            if par_sum >= rho:
                return x, y
            
    # se non restituisce dentro al ciclo for c'è qualcosa che non va
    # se mai succedesse, ritorno l'ultimo sito
    return Lx-1, Ly-1
