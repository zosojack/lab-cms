# KMCNumbaSubroutines.py
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
    Lx, Ly = height.shape
    
    for x in range(Lx):
        for y in range(Ly):
            first_neigh[x,y] = neigh_XY(x, y, height)
    
    return first_neigh

@njit
def _update_NN_deposition (first_neigh: np.ndarray,
                           k_diff: np.ndarray,                  # matrice dei rate di diffusione
                           current_k_diff_sum: float,           # somma prima dell'aggiornamento
                           current_k_diff_row_sums: np.ndarray, # somma per riga prima dell'aggiornamento
                           height: np.ndarray,
                           deposition_site: tuple[int, int],
                           rates_lookup: np.ndarray,
) -> np.ndarray:
    ''' 
    Aggiorna i primi vicini e la matrice dei rate di diffusione dopo un evento di deposizione;
    agisce solamente sui siti coinvolti nell'evento.
    '''
    x0, y0 = deposition_site
    Lx, Ly = height.shape
    
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
        n = int(first_neigh[x_dep, y_dep])
        # aggiorna rate di diffusione HACK: non ricalcolo ogni volta i rate
        current_k_diff_sum -= k_diff[x_dep, y_dep] # rimuovo il vecchio contributo
        current_k_diff_row_sums[x_dep] -= k_diff[x_dep, y_dep] # rimuovo il vecchio contributo
        
        k_diff[x_dep, y_dep] = rates_lookup[n]     # aggiorno con il nuovo rate
        
        current_k_diff_sum += k_diff[x_dep, y_dep] # aggiungo il nuovo contributo
        current_k_diff_row_sums[x_dep] += k_diff[x_dep, y_dep] # aggiungo il nuovo contributo
    
    return first_neigh, k_diff, current_k_diff_sum, current_k_diff_row_sums
    
@njit
def _update_NN_diffusion  (first_neigh: np.ndarray,
                           k_diff: np.ndarray,              # matrice dei rate di diffusione
                           current_k_diff_sum: float,       # somma prima dell'aggiornamento    
                           current_k_diff_row_sums: np.ndarray, # somma per riga prima dell'aggiornamento
                           height: np.ndarray,
                           diffusion_from: tuple[int, int],
                           diffusion_to: tuple[int, int],   
                           rates_lookup: np.ndarray,
) -> np.ndarray:
    ''' 
    Aggiorna i primi vicini e la matrice dei rate di diffusione dopo un evento di diffusione;
    agisce solamente sui siti coinvolti nell'evento.
    '''
    x_from, y_from = diffusion_from
    x_to, y_to = diffusion_to
    Lx, Ly = height.shape
    
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
        n = int(first_neigh[x_site, y_site])
        # aggiorna rate di diffusione HACK: non ricalcolo ogni volta i rate
        current_k_diff_sum -= k_diff[x_site, y_site] # rimuovo il vecchio contributo
        current_k_diff_row_sums[x_site] -= k_diff[x_site, y_site] # rimuovo il vecchio contributo
        
        k_diff[x_site, y_site] = rates_lookup[n]     # aggiorno con il nuovo rate
        
        current_k_diff_sum += k_diff[x_site, y_site] # aggiungo il nuovo contributo
        current_k_diff_row_sums[x_site] += k_diff[x_site, y_site] # aggiungo il nuovo contributo
    
    return first_neigh, k_diff, current_k_diff_sum, current_k_diff_row_sums

# selezione evento di diffusione
@njit
def find_move (rho: float,
               k_diff: np.ndarray,
               current_k_diff_row_sums: np.ndarray,
) -> tuple:
    ''' 
    Percorre la matrice dei rate di diffusione fino a raggiungere un valore >= rho;
    restituisce le coordinate (x,y) dell'atomo da muovere.
    Utilizza la somma per riga dei rate di diffusione per velocizzare la ricerca.
    '''
    Lx, Ly = k_diff.shape
    
    # 1. TROVA LA RIGA
    row_idx = 0
    current_sum = 0.0
    
    for x in range(Lx):
        r_sum = current_k_diff_row_sums[x]
        # Se sommando questa riga supero rho, l'evento è qui
        if current_sum + r_sum >= rho:
            row_idx = x
            break
        current_sum += r_sum
        
    # 2. CALCOLA IL RESIDUO
    # Rimuovo dal target le somme delle righe già scorse 
    rho -= current_sum 
    
    # 3. TROVA LA COLONNA
    col_sum = 0.0
    # Scorro in avanti
    for y in range(Ly):
        col_sum += k_diff[row_idx, y]
        if col_sum >= rho:
            return row_idx, y
            
    # Fallback di sicurezza
    return row_idx, Ly - 1