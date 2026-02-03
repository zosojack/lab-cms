# Kineticmontecarlo

[CMS Index](../../../README.md#cms-index) / `zosojack` / [Cms](../index.md#cms) / [Montecarlo](./index.md#montecarlo) / Kineticmontecarlo

> Auto-generated documentation for [zosojack.CMS.MonteCarlo.KineticMonteCarlo](../../../../zosojack/CMS/MonteCarlo/KineticMonteCarlo.py) module.

#### Attributes

- `k_B` - Costanti fisiche: 1 / 11603


- [Kineticmontecarlo](#kineticmontecarlo)
  - [KineticMonteCarlo](#kineticmontecarlo)
    - [KineticMonteCarlo()._initialize_diffusion_rates](#kineticmontecarlo()_initialize_diffusion_rates)
    - [KineticMonteCarlo().run](#kineticmontecarlo()run)
  - [KineticMonteCarloResult](#kineticmontecarloresult)

## KineticMonteCarlo

[Show source in KineticMonteCarlo.py:19](../../../../zosojack/CMS/MonteCarlo/KineticMonteCarlo.py#L19)

KineticMonteCarlo
=================
Classe per eseguire simulazioni kinetic Monte Carlo di crescita su reticolo.

Parameters
----------
L : tuple[int, int]
    Dimensioni del reticolo (Lx, Ly).
flux : float
    Flusso di deposizione in ML/s.
T : float
    Temperatura in Kelvin.
nu : float, optional
    Frequenza di tentativi di diffusione (default è 1.E+13 Hz).
J0 : float, optional
    Energia di interazione del primo vicino (default è 4 * J1 in eV).
J1 : float, optional
    Energia di interazione del secondo vicino (default è -0.345 eV).
seed : int, optional
    Seme per il generatore di numeri casuali (default è 123413432).

Attributes
----------
height : np.ndarray
    Matrice che rappresenta l'altezza della superficie in ogni sito del reticolo.
first_neigh : np.ndarray
    Matrice che contiene il numero di primi vicini per ciascun atomo sulla superficie.
k_diff : np.ndarray
    Matrice dei rate di diffusione per ciascun sito del reticolo.
k_diff_sum : float
    Somma totale dei rate di diffusione.

Methods
-------
run(end_time: float) -> KineticMonteCarloResult
    Esegue la simulazione fino al tempo specificato e restituisce i risultati
    raggruppati in un oggetto KineticMonteCarloResult.

#### Signature

```python
class KineticMonteCarlo:
    def __init__(
        self,
        L: tuple[int, int],
        flux: float,
        T: float,
        nu: float = 10000000000000.0,
        J0: float = 4 * -0.345,
        J1: float = -0.345,
        xyz_writer: Optional[XYZwriter] = None,
        seed: int = 123413432,
    ): ...
```

### KineticMonteCarlo()._initialize_diffusion_rates

[Show source in KineticMonteCarlo.py:102](../../../../zosojack/CMS/MonteCarlo/KineticMonteCarlo.py#L102)

Inizializza i primi vicini, i rate di diffusione e la loro somma

#### Signature

```python
def _initialize_diffusion_rates(self) -> None: ...
```

### KineticMonteCarlo().run

[Show source in KineticMonteCarlo.py:189](../../../../zosojack/CMS/MonteCarlo/KineticMonteCarlo.py#L189)

#### Signature

```python
def run(self, end_time: float) -> KineticMonteCarloResult: ...
```

#### See also

- [KineticMonteCarloResult](#kineticmontecarloresult)



## KineticMonteCarloResult

[Show source in KineticMonteCarlo.py:247](../../../../zosojack/CMS/MonteCarlo/KineticMonteCarlo.py#L247)

#### Attributes

- `n_deposition_events`: `int` - campi forniti in input

- `scarto_lunghezza_array`: `int` - campi calcolati in post-init: field(init=False)


KineticMonteCarloResult
=======================
Classe per gestire i risultati di una simulazione di kinetic Monte Carlo.

#### Signature

```python
class KineticMonteCarloResult: ...
```