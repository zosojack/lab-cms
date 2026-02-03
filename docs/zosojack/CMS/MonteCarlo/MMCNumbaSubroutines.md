# Mmcnumbasubroutines

[CMS Index](../../../README.md#cms-index) / `zosojack` / [Cms](../index.md#cms) / [Montecarlo](./index.md#montecarlo) / Mmcnumbasubroutines

> Auto-generated documentation for [zosojack.CMS.MonteCarlo.MMCNumbaSubroutines](../../../../zosojack/CMS/MonteCarlo/MMCNumbaSubroutines.py) module.

- [Mmcnumbasubroutines](#mmcnumbasubroutines)
  - [_compute_system_energy](#_compute_system_energy)
  - [_update_system_energy](#_update_system_energy)
  - [neigh_XYZ](#neigh_xyz)
  - [pbc_corr](#pbc_corr)

## _compute_system_energy

[Show source in MMCNumbaSubroutines.py:38](../../../../zosojack/CMS/MonteCarlo/MMCNumbaSubroutines.py#L38)

Ricalcola da zero l'energia del sistema

#### Signature

```python
@njit
def _compute_system_energy(height): ...
```



## _update_system_energy

[Show source in MMCNumbaSubroutines.py:58](../../../../zosojack/CMS/MonteCarlo/MMCNumbaSubroutines.py#L58)

Calcola la nuova energia totale aggiornando solo la differenza (Delta E).
Utilizza la matrice 'height' che è GIÀ stata modificata dalla trial move.

#### Signature

```python
@njit
def _update_system_energy(
    height: np.ndarray,
    current_energy: float,
    start_site: tuple[int, int],
    end_site: tuple[int, int],
) -> float: ...
```



## neigh_XYZ

[Show source in MMCNumbaSubroutines.py:16](../../../../zosojack/CMS/MonteCarlo/MMCNumbaSubroutines.py#L16)

Conta quanti primi vicini possiede l'atomo in posizione X,Y,Z

#### Signature

```python
@njit
def neigh_XYZ(x_i: int, y_i: int, z_i: int, height: np.ndarray) -> int: ...
```



## pbc_corr

[Show source in MMCNumbaSubroutines.py:5](../../../../zosojack/CMS/MonteCarlo/MMCNumbaSubroutines.py#L5)

Correzione condizione periodica al contorno

#### Signature

```python
@njit
def pbc_corr(xi: int, Li: int) -> int: ...
```