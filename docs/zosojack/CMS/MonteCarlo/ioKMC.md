# Iokmc

[CMS Index](../../../README.md#cms-index) / `zosojack` / [Cms](../index.md#cms) / [Montecarlo](./index.md#montecarlo) / Iokmc

> Auto-generated documentation for [zosojack.CMS.MonteCarlo.ioKMC](../../../../zosojack/CMS/MonteCarlo/ioKMC.py) module.

- [Iokmc](#iokmc)
  - [XYZwriter](#xyzwriter)
    - [XYZwriter().write_frame](#xyzwriter()write_frame)
  - [_convert_pile_to_xyz](#_convert_pile_to_xyz)

## XYZwriter

[Show source in ioKMC.py:36](../../../../zosojack/CMS/MonteCarlo/ioKMC.py#L36)

XYZwriter
=========

Classe per registrare le traiettorie atomiche in formato XYZ.

#### Signature

```python
class XYZwriter:
    def __init__(
        self, output_folder: str, max_deposition_steps: int = 1000, max_files: int = 500
    ) -> None: ...
```

### XYZwriter().write_frame

[Show source in ioKMC.py:60](../../../../zosojack/CMS/MonteCarlo/ioKMC.py#L60)

Costruisce e registra le posizioni degli atomi depositati sulla superficie.
L'azione è eseguita un massimo di max_files volte.

Parameters
----------
time : float
    Tempo corrente della simulazione.
deposition_step : int
    Passo temporale corrente della simulazione.
height : np.ndarray
    Altezze delle pile di atomi sulla superficie.

#### Signature

```python
def write_frame(self, time: float, deposition_step: int, height: np.ndarray) -> None: ...
```



## _convert_pile_to_xyz

[Show source in ioKMC.py:6](../../../../zosojack/CMS/MonteCarlo/ioKMC.py#L6)

Converte una matrice di altezze (pile di atomi) in una matrice di posizioni (N_atoms x 3).
Ogni riga della matrice risultante rappresenta le coordinate (x, y, z) di un atomo.

#### Signature

```python
@njit
def _convert_pile_to_xyz(height: np.ndarray) -> np.ndarray: ...
```