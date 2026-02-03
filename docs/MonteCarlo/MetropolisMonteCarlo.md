# Metropolismontecarlo

[CMS Index](../README.md#cms-index) / [Montecarlo](./index.md#montecarlo) / Metropolismontecarlo

> Auto-generated documentation for [MonteCarlo.MetropolisMonteCarlo](../../zosojack/CMS/MonteCarlo/MetropolisMonteCarlo.py) module.

- [Metropolismontecarlo](#metropolismontecarlo)
  - [MetropolisMonteCarlo](#metropolismontecarlo)

## MetropolisMonteCarlo

[Show source in MetropolisMonteCarlo.py:5](../../zosojack/CMS/MonteCarlo/MetropolisMonteCarlo.py#L5)

MetropolisMonteCarlo
====================
Classe per eseguire simulazioni Monte Carlo di tipo Metropolis su reticolo.

Parameters
----------
L : tuple[int, int]
    Dimensioni del reticolo (Lx, Ly).
T : float
    Temperatura in Kelvin.
J : float
    Energia di interazione tra spin (default è 1.0 eV).
seed : int, optional
    Seme per il generatore di numeri casuali (default è 123413432).

Attributes
----------
spins : np.ndarray
    Matrice che rappresenta gli spin in ogni sito del reticolo.

Methods
-------
run(steps: int) -> MetropolisMonteCarloResult
    Esegue la simulazione per un numero specificato di passi e restituisce i risultati
    raggruppati in un oggetto MetropolisMonteCarloResult.

#### Signature

```python
class MetropolisMonteCarlo:
    def __init__(
        self, L: tuple[int, int], T: float, J: float = 1.0, seed: int = 123413432
    ): ...
```