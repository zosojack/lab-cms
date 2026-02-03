# Iommc

[CMS Index](../../../README.md#cms-index) / `zosojack` / [Cms](../index.md#cms) / [Montecarlo](./index.md#montecarlo) / Iommc

> Auto-generated documentation for [zosojack.CMS.MonteCarlo.ioMMC](../../../../zosojack/CMS/MonteCarlo/ioMMC.py) module.

- [Iommc](#iommc)
  - [XYZwriter](#xyzwriter)
    - [XYZwriter().write_frame](#xyzwriter()write_frame)
  - [_convert_pile_to_xyz](#_convert_pile_to_xyz)

## XYZwriter

[Show source in ioMMC.py:41](../../../../zosojack/CMS/MonteCarlo/ioMMC.py#L41)

XYZwriter
=========

Classe per registrare le traiettorie atomiche in formato XYZ.

#### Signature

```python
class XYZwriter:
    def __init__(self, output_folder: str, max_files: int = 500) -> None: ...
```

### XYZwriter().write_frame

[Show source in ioMMC.py:61](../../../../zosojack/CMS/MonteCarlo/ioMMC.py#L61)

Costruisce e registra le posizioni degli atomi depositati sulla superficie.
L'azione è eseguita un massimo di max_files volte.

Parameters
----------
current_accepted_step : int
    Numero di volte che è stata accettata una mossa della simulazione.
height : np.ndarray
    Altezze delle pile di atomi sulla superficie.

#### Signature

```python
def write_frame(self, current_accepted_step: int, height: np.ndarray) -> None: ...
```



## _convert_pile_to_xyz

[Show source in ioMMC.py:12](../../../../zosojack/CMS/MonteCarlo/ioMMC.py#L12)

Converte una matrice di altezze (pile di atomi) in una matrice di posizioni (N_atoms x 3).
Ogni riga della matrice risultante rappresenta le coordinate (x, y, z) di un atomo.

#### Signature

```python
@njit
def _convert_pile_to_xyz(height: np.ndarray) -> np.ndarray: ...
```