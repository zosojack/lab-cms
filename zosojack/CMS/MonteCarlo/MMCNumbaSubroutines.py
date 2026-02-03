# MMCNumbaSubroutines.py
import numpy as np
from numba import njit

@njit
def pbc_corr (xi: int, 
              Li: int
) -> int:
    ''' Correzione condizione periodica al contorno '''
    
    return xi % Li
    # Alternative implementation:
    xi -= Li * np.floor(xi/Li)
    return int(xi)

@njit
def neigh_XYZ (x_i: int,
               y_i: int,
               z_i: int,
               height: np.ndarray # potrei passare solo una slice della matrice, invece che tutta
) -> int:
    ''' Conta quanti primi vicini possiede l'atomo in posizione X,Y,Z '''
    
    Lx, Ly = height.shape[0], height.shape[1]
    n = 0
    
    # Vicino 1 (Up)
    if height[x_i, pbc_corr(y_i+1, Ly)] >= z_i: n += 1
    # Vicino 2 (Down)
    if height[x_i, pbc_corr(y_i-1, Ly)] >= z_i: n += 1
    # Vicino 3 (Right)
    if height[pbc_corr(x_i+1, Lx), y_i] >= z_i: n += 1
    # Vicino 4 (Left)
    if height[pbc_corr(x_i-1, Lx), y_i] >= z_i: n += 1
            
    return int(n)

@njit
def _compute_system_energy(height):
    ''' Ricalcola da zero l'energia del sistema '''
    en=0
    Lx, Ly = height.shape
    J0 = 4*(-0.345) # default: J0 = 4 * J1 [eV]
    J1 = -0.345  # energia di legame orizzontale      
    J1_2 = J1 / 2.0
    for x in range(Lx):
        for y in range(Ly):
            # QUI POTREBBE ESSERE INTELLIGENTE PASSARE UNA SLICE DELLA MATRICE height
            for z in range(1,height[x,y]+1):
                n1=neigh_XYZ(x,y,z,height)
                en = en + J0 + J1_2*n1
    return en




# TODO: controllare questa funzione
@njit
def _update_system_energy(height: np.ndarray, 
                          current_energy: float,
                          start_site: tuple[int, int],
                          end_site: tuple[int, int]) -> float:
    ''' 
    Calcola la nuova energia totale aggiornando solo la differenza (Delta E).
    Utilizza la matrice 'height' che è GIÀ stata modificata dalla trial move.
    '''
    
    # Estraiamo le coordinate (ignoriamo le z passate nella tupla perché 
    # ricalcoliamo le altezze reali dalla matrice per sicurezza)
    x_from, y_from = start_site
    x_to, y_to = end_site
    
    # Costanti dell'esercizio (DEVONO essere identiche a compute_system_energy)
    J1 = -0.345          # Energia laterale
    # J0 si elide nel calcolo del Delta E, quindi non serve esplicitarlo
    
    # 1. Calcoliamo i legami nel sito di ARRIVO (dove l'atomo è ORA presente)
    z_new = height[x_to, y_to]
    n_new = neigh_XYZ(x_to, y_to, z_new, height)
    
    # 2. "Flip-Flop": Ripristiniamo temporaneamente la matrice allo stato precedente
    #    per vedere quanti legami aveva l'atomo nel sito di partenza.
    height[x_to, y_to] -= 1
    height[x_from, y_from] += 1
    
    # 3. Calcoliamo i legami nel sito di PARTENZA (dove l'atomo ERA presente)
    z_old = height[x_from, y_from]
    n_old = neigh_XYZ(x_from, y_from, z_old, height)
    
    # 4. Ripristiniamo la modifica (rimettiamo l'atomo nel sito di arrivo come era all'inizio della funzione)
    height[x_from, y_from] -= 1
    height[x_to, y_to] += 1
    
    # 5. Calcolo Delta E
    # Delta = (Legami Nuovi - Legami Vecchi) * J1
    # Usiamo J1 intero (non mezzi) perché stiamo spostando un intero atomo e i suoi legami.
    delta_E = (n_new - n_old) * J1
    
    return current_energy + delta_E