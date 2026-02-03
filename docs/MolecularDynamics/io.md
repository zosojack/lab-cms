# Io

[CMS Index](../README.md#cms-index) / [Moleculardynamics](./index.md#moleculardynamics) / Io

> Auto-generated documentation for [MolecularDynamics.io](../../zosojack/CMS/MolecularDynamics/io.py) module.

- [Io](#io)
  - [AtomTracker](#atomtracker)
    - [AtomTracker().record_position](#atomtracker()record_position)
    - [AtomTracker().set_pcb_option](#atomtracker()set_pcb_option)
  - [XYZwriter](#xyzwriter)
    - [XYZwriter().write_frame](#xyzwriter()write_frame)

## AtomTracker

[Show source in io.py:9](../../zosojack/CMS/MolecularDynamics/io.py#L9)

AtomTracker
===========
Classe che traccia le posizioni di uno specifico atomo durante la simulazione.

Attributes
----------
index : int
    Indice dell'atomo da tracciare.
output_file : str
    Percorso del file di output dove salvare le posizioni.
pbc_option : str
    Modalità per le condizioni al contorno: 'periodic' oppure 'unbounded' (default).

Methods
-------
record_position(step: int, crystal: CrystalStructure) -> None
    Registra la posizione dell'atomo al passo specificato.

#### Signature

```python
class AtomTracker:
    def __init__(
        self, index: int, output_file: str, pbc_option: str = "unbounded"
    ) -> None: ...
```

### AtomTracker().record_position

[Show source in io.py:58](../../zosojack/CMS/MolecularDynamics/io.py#L58)

Registra le posizioni degli atomi tracciati nel file di output.

Parameters
----------
step : int
    Passo temporale corrente della simulazione.
crystal : CrystalStructure
    Struttura cristallina con posizioni e dimensioni della cella per le PBC.

#### Signature

```python
def record_position(self, step: int, crystal: CrystalStructure) -> None: ...
```

### AtomTracker().set_pcb_option

[Show source in io.py:45](../../zosojack/CMS/MolecularDynamics/io.py#L45)

Imposta l'opzione per il trattamento delle condizioni al contorno periodiche.

Parameters
----------
option : str
    'periodic' per condizioni al contorno periodiche, 'unbounded' altrimenti.

#### Signature

```python
def set_pcb_option(self, option: str) -> None: ...
```



## XYZwriter

[Show source in io.py:78](../../zosojack/CMS/MolecularDynamics/io.py#L78)

XYZwriter
==========
Classe che registra le traiettorie atomiche in formato .xyz.

Attributes
----------
output_folder : str
    Cartella di output dove salvare le traiettorie.
dt : float
    Intervallo di tempo tra i frame registrati.
dump_interval : int
    Passi tra i frame salvati (default: 200).

Methods
-------
write_frame(step: int, positions: np.ndarray) -> None
    Registra le posizioni di tutti gli atomi nel file di output.

#### Signature

```python
class XYZwriter:
    def __init__(
        self, output_folder: str, dt: float, dump_interval: int = 200
    ) -> None: ...
```

### XYZwriter().write_frame

[Show source in io.py:108](../../zosojack/CMS/MolecularDynamics/io.py#L108)

Registra le posizioni di tutti gli atomi nel file di output.
L'azione è eseguita soltanto ogni dump_interval steps.

Parameters
----------
positions : np.ndarray
    Posizioni attuali di tutti gli atomi.
step : int
    Passo temporale corrente della simulazione.

#### Signature

```python
def write_frame(self, step: int, positions: np.ndarray) -> None: ...
```