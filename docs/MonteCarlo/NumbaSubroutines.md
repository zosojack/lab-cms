# Numbasubroutines

[CMS Index](../README.md#cms-index) / [Montecarlo](./index.md#montecarlo) / Numbasubroutines

> Auto-generated documentation for [MonteCarlo.NumbaSubroutines](../../zosojack/CMS/MonteCarlo/NumbaSubroutines.py) module.

- [Numbasubroutines](#numbasubroutines)
  - [_update_NN_deposition](#_update_nn_deposition)
  - [_update_NN_diffusion](#_update_nn_diffusion)
  - [count_NN](#count_nn)
  - [find_move](#find_move)
  - [neigh_XY](#neigh_xy)
  - [pbc_corr](#pbc_corr)

## _update_NN_deposition

[Show source in NumbaSubroutines.py:50](../../zosojack/CMS/MonteCarlo/NumbaSubroutines.py#L50)

Aggiorna i primi vicini e la matrice dei rate di diffusione dopo un evento di deposizione;
agisce solamente sui siti coinvolti nell'evento.

#### Signature

```python
@njit
def _update_NN_deposition(
    first_neigh: np.ndarray,
    k_diff: np.ndarray,
    current_k_diff_sum: float,
    height: np.ndarray,
    deposition_site: tuple[int, int],
    params: tuple,
) -> np.ndarray: ...
```



## _update_NN_diffusion

[Show source in NumbaSubroutines.py:88](../../zosojack/CMS/MonteCarlo/NumbaSubroutines.py#L88)

Aggiorna i primi vicini e la matrice dei rate di diffusione dopo un evento di diffusione;
agisce solamente sui siti coinvolti nell'evento.

#### Signature

```python
@njit
def _update_NN_diffusion(
    first_neigh: np.ndarray,
    k_diff: np.ndarray,
    current_k_diff_sum: float,
    height: np.ndarray,
    diffusion_from: tuple[int, int],
    diffusion_to: tuple[int, int],
    params: tuple,
) -> np.ndarray: ...
```



## count_NN

[Show source in NumbaSubroutines.py:37](../../zosojack/CMS/MonteCarlo/NumbaSubroutines.py#L37)

Produce una matrice di interi che contiene il numero di primi vicini per ciascun atomo sulla superficie

#### Signature

```python
@njit
def count_NN(height: np.ndarray) -> np.ndarray: ...
```



## find_move

[Show source in NumbaSubroutines.py:134](../../zosojack/CMS/MonteCarlo/NumbaSubroutines.py#L134)

Percorre la matrice dei rate di diffusione fino a raggiungere un valore >= rho;
restituisce le coordinate (x,y) dell'atomo da muovere.

#### Signature

```python
@njit
def find_move(rho: float, k_diff: np.ndarray) -> tuple: ...
```



## neigh_XY

[Show source in NumbaSubroutines.py:15](../../zosojack/CMS/MonteCarlo/NumbaSubroutines.py#L15)

Conta quanti primi vicini possiede l'atomo sulla superficie in posizione X,Y

#### Signature

```python
@njit
def neigh_XY(x_i: float, y_i: float, height: np.ndarray) -> int: ...
```



## pbc_corr

[Show source in NumbaSubroutines.py:6](../../zosojack/CMS/MonteCarlo/NumbaSubroutines.py#L6)

Correzione condizione periodica al contorno

#### Signature

```python
@njit
def pbc_corr(xi: float, Li: float) -> int: ...
```