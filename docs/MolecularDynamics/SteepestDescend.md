# Steepestdescend

[CMS Index](../README.md#cms-index) / [Moleculardynamics](./index.md#moleculardynamics) / Steepestdescend

> Auto-generated documentation for [MolecularDynamics.SteepestDescend](../../zosojack/CMS/MolecularDynamics/SteepestDescend.py) module.

- [Steepestdescend](#steepestdescend)
  - [SteepestDescend](#steepestdescend)
    - [SteepestDescend().minimize_energy](#steepestdescend()minimize_energy)

## SteepestDescend

[Show source in SteepestDescend.py:9](../../zosojack/CMS/MolecularDynamics/SteepestDescend.py#L9)

SteepestDescent
===============
Classe che implementa la minimizzazione dell'energia tramite il metodo del gradiente discendente.

Attributes
----------
crystal : CrystalStructure
    Struttura cristallina da minimizzare.

Methods
-------
minimize_energy(max_steps=1000, F_tol=1e-5, C_steep=0.01, pol_junction=False) -> Tuple[np.ndarray, np.ndarray]
    Esegue la minimizzazione dell'energia e restituisce l'energia potenziale e le forze massime ad ogni passo.

#### Signature

```python
class SteepestDescend:
    def __init__(self, crystal: CrystalStructure): ...
```

### SteepestDescend().minimize_energy

[Show source in SteepestDescend.py:31](../../zosojack/CMS/MolecularDynamics/SteepestDescend.py#L31)

Esegue la minimizzazione dell'energia e restituisce l'energia potenziale e le forze massime ad ogni passo.

Parameters
----------
max_steps : int, optional
    Numero massimo di passi per la minimizzazione (default è 1000).
F_tol : float, optional
    Tolleranza per la forza massima per la convergenza (default è 1e-5).
C_steep : float, optional
    Coefficiente di discesa per l'aggiornamento delle posizioni (default è 0.01).
pol_junction : bool, optional
    Se True, utilizza il giunto polinomiale per il potenziale (default è False).

#### Signature

```python
def minimize_energy(
    self, max_steps=1000, F_tol=1e-05, C_steep=0.01, pol_junction=False
): ...
```