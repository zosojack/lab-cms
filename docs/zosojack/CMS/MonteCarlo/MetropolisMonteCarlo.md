# Metropolismontecarlo

[CMS Index](../../../README.md#cms-index) / `zosojack` / [Cms](../index.md#cms) / [Montecarlo](./index.md#montecarlo) / Metropolismontecarlo

> Auto-generated documentation for [zosojack.CMS.MonteCarlo.MetropolisMonteCarlo](../../../../zosojack/CMS/MonteCarlo/MetropolisMonteCarlo.py) module.

#### Attributes

- `k_B` - Costanti fisiche: 1 / 11603


- [Metropolismontecarlo](#metropolismontecarlo)
  - [MetropolisMonteCarlo](#metropolismontecarlo)
    - [MetropolisMonteCarlo()._select_end_site](#metropolismontecarlo()_select_end_site)
    - [MetropolisMonteCarlo()._select_start_site](#metropolismontecarlo()_select_start_site)
    - [MetropolisMonteCarlo().run](#metropolismontecarlo()run)
  - [MetropolisMonteCarloResult](#metropolismontecarloresult)

## MetropolisMonteCarlo

[Show source in MetropolisMonteCarlo.py:14](../../../../zosojack/CMS/MonteCarlo/MetropolisMonteCarlo.py#L14)

MetropolisMonteCarlo
====================
Classe per eseguire simulazioni Monte Carlo di tipo Metropolis su reticolo.
Può individuare la configurazione di equilibrio di un sistema di atomi su un reticolo 2D,
non prevederne le traiettorie. Non contiene variabili dipendenti dal tempo.

Parameters
----------
L : tuple[int, int]
    Dimensioni del reticolo (Lx, Ly).
N_atoms : int
    Numero totale di atomi nel reticolo.
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
        self,
        L: tuple[int, int],
        N_atoms: int,
        T: float,
        J: float = 1.0,
        seed: int = 123413432,
        xyz_writer: Optional[XYZwriter] = None,
    ): ...
```

### MetropolisMonteCarlo()._select_end_site

[Show source in MetropolisMonteCarlo.py:89](../../../../zosojack/CMS/MonteCarlo/MetropolisMonteCarlo.py#L89)

Seleziona un sito casuale tra tutti i siti del reticolo

#### Signature

```python
def _select_end_site(self) -> tuple[int, int]: ...
```

### MetropolisMonteCarlo()._select_start_site

[Show source in MetropolisMonteCarlo.py:84](../../../../zosojack/CMS/MonteCarlo/MetropolisMonteCarlo.py#L84)

Seleziona un sito casuale tra quelli occupati

#### Signature

```python
def _select_start_site(self) -> tuple[int, int]: ...
```

### MetropolisMonteCarlo().run

[Show source in MetropolisMonteCarlo.py:97](../../../../zosojack/CMS/MonteCarlo/MetropolisMonteCarlo.py#L97)

#### Signature

```python
def run(
    self, N_steps: int = 100000, thermalization_steps: int = 0
) -> MetropolisMonteCarloResult: ...
```

#### See also

- [MetropolisMonteCarloResult](#metropolismontecarloresult)



## MetropolisMonteCarloResult

[Show source in MetropolisMonteCarlo.py:170](../../../../zosojack/CMS/MonteCarlo/MetropolisMonteCarlo.py#L170)

MetropolisMonteCarloResult
==========================
Classe per memorizzare i risultati di una simulazione Metropolis Monte Carlo.

Parameters
----------
energies : np.ndarray
    Lista delle energie totali ad ogni passo.
acceptance_ratio : float
    Rapporto di accettazione degli spostamenti proposti.

#### Signature

```python
class MetropolisMonteCarloResult: ...
```