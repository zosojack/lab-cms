# Kmcnumbasubroutines

[CMS Index](../../../README.md#cms-index) / `zosojack` / [Cms](../index.md#cms) / [Montecarlo](./index.md#montecarlo) / Kmcnumbasubroutines

> Auto-generated documentation for [zosojack.CMS.MonteCarlo.KMCNumbaSubroutines](../../../../zosojack/CMS/MonteCarlo/KMCNumbaSubroutines.py) module.

- [Kmcnumbasubroutines](#kmcnumbasubroutines)
  - [_update_NN_deposition](#_update_nn_deposition)
  - [_update_NN_diffusion](#_update_nn_diffusion)
  - [count_NN](#count_nn)
  - [find_move](#find_move)
  - [neigh_XY](#neigh_xy)
  - [pbc_corr](#pbc_corr)

## _update_NN_deposition

[Show source in KMCNumbaSubroutines.py:50](../../../../zosojack/CMS/MonteCarlo/KMCNumbaSubroutines.py#L50)

Aggiorna i primi vicini e la matrice dei rate di diffusione dopo un evento di deposizione;
agisce solamente sui siti coinvolti nell'evento.

#### Signature

```python
@njit
def _update_NN_deposition(
    first_neigh: np.ndarray,
    k_diff: np.ndarray,
    current_k_diff_sum: float,
    current_k_diff_row_sums: np.ndarray,
    height: np.ndarray,
    deposition_site: tuple[int, int],
    rates_lookup: np.ndarray,
) -> np.ndarray: ...
```



## _update_NN_diffusion

[Show source in KMCNumbaSubroutines.py:91](../../../../zosojack/CMS/MonteCarlo/KMCNumbaSubroutines.py#L91)

Aggiorna i primi vicini e la matrice dei rate di diffusione dopo un evento di diffusione;
agisce solamente sui siti coinvolti nell'evento.

#### Signature

```python
@njit
def _update_NN_diffusion(
    first_neigh: np.ndarray,
    k_diff: np.ndarray,
    current_k_diff_sum: float,
    current_k_diff_row_sums: np.ndarray,
    height: np.ndarray,
    diffusion_from: tuple[int, int],
    diffusion_to: tuple[int, int],
    rates_lookup: np.ndarray,
) -> np.ndarray: ...
```



## count_NN

[Show source in KMCNumbaSubroutines.py:37](../../../../zosojack/CMS/MonteCarlo/KMCNumbaSubroutines.py#L37)

Produce una matrice di interi che contiene il numero di primi vicini per ciascun atomo sulla superficie

#### Signature

```python
@njit
def count_NN(height: np.ndarray) -> np.ndarray: ...
```



## find_move

[Show source in KMCNumbaSubroutines.py:140](../../../../zosojack/CMS/MonteCarlo/KMCNumbaSubroutines.py#L140)

Percorre la matrice dei rate di diffusione fino a raggiungere un valore >= rho;
restituisce le coordinate (x,y) dell'atomo da muovere.
Utilizza la somma per riga dei rate di diffusione per velocizzare la ricerca.

#### Signature

```python
@njit
def find_move(
    rho: float, k_diff: np.ndarray, current_k_diff_row_sums: np.ndarray
) -> tuple: ...
```



## neigh_XY

[Show source in KMCNumbaSubroutines.py:15](../../../../zosojack/CMS/MonteCarlo/KMCNumbaSubroutines.py#L15)

Conta quanti primi vicini possiede l'atomo sulla superficie in posizione X,Y

#### Signature

```python
@njit
def neigh_XY(x_i: float, y_i: float, height: np.ndarray) -> int: ...
```



## pbc_corr

[Show source in KMCNumbaSubroutines.py:6](../../../../zosojack/CMS/MonteCarlo/KMCNumbaSubroutines.py#L6)

Correzione condizione periodica al contorno

#### Signature

```python
@njit
def pbc_corr(xi: float, Li: float) -> int: ...
```