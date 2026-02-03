# Crystaldynamics

[CMS Index](../README.md#cms-index) / [Moleculardynamics](./index.md#moleculardynamics) / Crystaldynamics

> Auto-generated documentation for [MolecularDynamics.CrystalDynamics](../../zosojack/CMS/MolecularDynamics/CrystalDynamics.py) module.

#### Attributes

- `k_B` - Costanti fisiche: 1 / 11603


- [Crystaldynamics](#crystaldynamics)
  - [CrystalDynamics](#crystaldynamics)
    - [CrystalDynamics()._output_positions](#crystaldynamics()_output_positions)
    - [CrystalDynamics()._output_state](#crystaldynamics()_output_state)
    - [CrystalDynamics()._random_velocities](#crystaldynamics()_random_velocities)
    - [CrystalDynamics()._temperature](#crystaldynamics()_temperature)
    - [CrystalDynamics()._update_forces](#crystaldynamics()_update_forces)
    - [CrystalDynamics()._update_neighbours_distances](#crystaldynamics()_update_neighbours_distances)
    - [CrystalDynamics()._update_positions](#crystaldynamics()_update_positions)
    - [CrystalDynamics()._update_velocities](#crystaldynamics()_update_velocities)
    - [CrystalDynamics().atom_tracker](#crystaldynamics()atom_tracker)
    - [CrystalDynamics().atom_tracker](#crystaldynamics()atom_tracker-1)
    - [CrystalDynamics().run_dynamics](#crystaldynamics()run_dynamics)
    - [CrystalDynamics().set_seed](#crystaldynamics()set_seed)
    - [CrystalDynamics().xyz_writer](#crystaldynamics()xyz_writer)
    - [CrystalDynamics().xyz_writer](#crystaldynamics()xyz_writer-1)
  - [CrystalDynamicsResult](#crystaldynamicsresult)
    - [CrystalDynamicsResult().__post_init__](#crystaldynamicsresult()__post_init__)
    - [CrystalDynamicsResult().summary](#crystaldynamicsresult()summary)

## CrystalDynamics

[Show source in CrystalDynamics.py:19](../../zosojack/CMS/MolecularDynamics/CrystalDynamics.py#L19)

CrystalDynamics
===============
Classe per eseguire la dinamica molecolare classica su un sistema cristallino.

Attributes
----------
crystal : CrystalStructure
 Oggetto CrystalStructure che rappresenta la struttura cristallina.
atomic_mass : float
 Massa atomica in unità di massa atomica (default: 108 u).
dt : float
 Passo temporale della simulazione in secondi (default: 1e-15 s).
temp_ini : float
 Temperatura iniziale in Kelvin (default: 20 K).
atom_tracker : AtomTracker | list[AtomTracker] | None
 Strumento(i) per tracciare la posizione di specifici atomi (opzionale).
xyz_writer : XYZwriter | None
 Strumento per salvare le posizioni in formato XYZ (opzionale).

Methods
-------
run_dynamics(n_steps: int, t_th: float = 0, rescale_velocity: bool = False, debug: bool = False)      -> CrystalDynamicsResult
 Esegue la dinamica molecolare per un numero specificato di passi.

#### Signature

```python
class CrystalDynamics:
    def __init__(
        self,
        crystal: CrystalStructure,
        atomic_mass: float = 108,
        dt: float = 1e-15,
        temp_ini: float = 20,
        atom_tracker: Optional[AtomTracker | list[AtomTracker]] = None,
        xyz_writer: Optional[XYZwriter] = None,
    ): ...
```

### CrystalDynamics()._output_positions

[Show source in CrystalDynamics.py:189](../../zosojack/CMS/MolecularDynamics/CrystalDynamics.py#L189)

Salva le posizioni istantanee (stesso formato dell'originale).

#### Signature

```python
def _output_positions(self, foldername, step, n_steps) -> None: ...
```

### CrystalDynamics()._output_state

[Show source in CrystalDynamics.py:182](../../zosojack/CMS/MolecularDynamics/CrystalDynamics.py#L182)

Salva lo stato della simulazione.

#### Signature

```python
def _output_state(self, filename, step, E_tot, E_pot, E_kin, temp) -> None: ...
```

### CrystalDynamics()._random_velocities

[Show source in CrystalDynamics.py:108](../../zosojack/CMS/MolecularDynamics/CrystalDynamics.py#L108)

Inizializza le velocità casuali (stesso algoritmo dell'originale).

#### Signature

```python
def _random_velocities(self) -> None: ...
```

### CrystalDynamics()._temperature

[Show source in CrystalDynamics.py:173](../../zosojack/CMS/MolecularDynamics/CrystalDynamics.py#L173)

Restituisce la temperatura attuale del sistema.

#### Signature

```python
def _temperature(self) -> float: ...
```

### CrystalDynamics()._update_forces

[Show source in CrystalDynamics.py:148](../../zosojack/CMS/MolecularDynamics/CrystalDynamics.py#L148)

Calcola le forze usando la versione numpy del potenziale.

#### Signature

```python
def _update_forces(self, poly7=None) -> None: ...
```

### CrystalDynamics()._update_neighbours_distances

[Show source in CrystalDynamics.py:132](../../zosojack/CMS/MolecularDynamics/CrystalDynamics.py#L132)

Aggiorna le distanze tra gli atomi.
Se le posizioni non sono cambiate più di un valore soglia, sono ricalcolate soltanto
le distanze tra i vicini già noti con only_neighbours_distance().
Se le posizioni sono cambiate troppo, il flag neighbours_computed è posto False e
viene rieseguito find_neighbours().

#### Signature

```python
def _update_neighbours_distances(self) -> None: ...
```

### CrystalDynamics()._update_positions

[Show source in CrystalDynamics.py:125](../../zosojack/CMS/MolecularDynamics/CrystalDynamics.py#L125)

Aggiorna le posizioni con Velocity-Verlet (stessa formula).

#### Signature

```python
def _update_positions(self) -> None: ...
```

### CrystalDynamics()._update_velocities

[Show source in CrystalDynamics.py:153](../../zosojack/CMS/MolecularDynamics/CrystalDynamics.py#L153)

Aggiorna le velocità.
Se rescale_velocity è True, ricalcola la velocità per mantenere la temperatura costante.

#### Signature

```python
def _update_velocities(self, rescale_velocity: bool = False) -> None: ...
```

### CrystalDynamics().atom_tracker

[Show source in CrystalDynamics.py:92](../../zosojack/CMS/MolecularDynamics/CrystalDynamics.py#L92)

#### Signature

```python
@property
def atom_tracker(self): ...
```

### CrystalDynamics().atom_tracker

[Show source in CrystalDynamics.py:96](../../zosojack/CMS/MolecularDynamics/CrystalDynamics.py#L96)

#### Signature

```python
@atom_tracker.setter
def atom_tracker(self, tracker: AtomTracker | list[AtomTracker]) -> None: ...
```

### CrystalDynamics().run_dynamics

[Show source in CrystalDynamics.py:202](../../zosojack/CMS/MolecularDynamics/CrystalDynamics.py#L202)

Esegue la dinamica molecolare per `n_steps` step.

Parameters
----------
n_steps : int
 Numero totale di passi della simulazione.
t_th : float, optional
 Tempo di termalizzazione in secondi (default: 0 s).
rescale_velocity : bool, optional
 Se True, ricalcola la velocità per mantenere la temperatura costante (default: False).
debug : bool, optional
 Se True, stampa informazioni di debug ad ogni passo (default: False).
track_last : bool, optional
 Se True, traccia l'ultimo atomo (default: False).
n_print : int | None, optional
 Numero di passi tra le stampe di debug (default: None).
output : bool, optional
 Se True, salva lo stato della simulazione ad ogni passo (default: False).

Returns
-------
CrystalDynamicsResult
 Oggetto contenente i risultati della simulazione.

#### Signature

```python
def run_dynamics(
    self,
    n_steps: float,
    t_th: float = 0,
    rescale_velocity: bool = False,
    debug: bool = False,
    track_last=False,
    n_print=None,
    output=False,
) -> CrystalDynamicsResult: ...
```

#### See also

- [CrystalDynamicsResult](#crystaldynamicsresult)

### CrystalDynamics().set_seed

[Show source in CrystalDynamics.py:104](../../zosojack/CMS/MolecularDynamics/CrystalDynamics.py#L104)

Imposta il seed del generatore pseudo-casuale.

#### Signature

```python
def set_seed(self, myseed) -> None: ...
```

### CrystalDynamics().xyz_writer

[Show source in CrystalDynamics.py:80](../../zosojack/CMS/MolecularDynamics/CrystalDynamics.py#L80)

#### Signature

```python
@property
def xyz_writer(self): ...
```

### CrystalDynamics().xyz_writer

[Show source in CrystalDynamics.py:84](../../zosojack/CMS/MolecularDynamics/CrystalDynamics.py#L84)

#### Signature

```python
@xyz_writer.setter
def xyz_writer(self, writer: XYZwriter) -> None: ...
```



## CrystalDynamicsResult

[Show source in CrystalDynamics.py:327](../../zosojack/CMS/MolecularDynamics/CrystalDynamics.py#L327)

#### Attributes

- `time_step`: `float` - --- 1. CAMPI OBBLIGATORI (vanno passati quando crei l'oggetto) ---

- `trajectory_folder_path`: `Optional[str]` - --- 2. CAMPI OPZIONALI (hanno un valore di default) ---: None

- `metadata`: `Dict` - Per liste/dict vuoti di default, si usa field(default_factory=...): field(default_factory=dict)

- `mean_temp`: `float` - --- 3. CAMPI CALCOLATI (non li passi, li calcola lui) ---
  Usiamo init=False per dire "non chiederlo nel costruttore": field(init=False)


 CrystalDynamicsResult
=====================
 Classe per gestire i risultati di una simulazione di dinamica molecolare.

Attributes
----------
time_step : float
 Il passo temporale della simulazione in secondi.
num_steps : int
 Il numero totale di passi della simulazione.
steps_th : int
 Il numero di passi di termalizzazione.
energies : Dict[str, np.ndarray]
 Un dizionario contenente array di energie ('total', 'kinetic', 'potential').
temperatures : np.ndarray
 Un array contenente le temperature ad ogni passo.
trajectory_folder_path : Optional[str]
 Il percorso della cartella contenente i file di traiettoria XYZ (se salvati).
adatom_file_path : Optional[str]
 Una lista di percorsi dei file di tracking degli adatom (se salvati).
mean_temp : float
 La temperatura media calcolata dopo la simulazione.
final_temp : float
 La temperatura finale della simulazione.
mean_E_tot : float
 L'energia totale media calcolata dopo la simulazione.
std_E_tot : float
 La deviazione standard dell'energia totale calcolata dopo la simulazione.

Methods
-------
summary() -> None
 Stampa un riepilogo dei risultati della simulazione.

#### Signature

```python
class CrystalDynamicsResult: ...
```

### CrystalDynamicsResult().__post_init__

[Show source in CrystalDynamics.py:385](../../zosojack/CMS/MolecularDynamics/CrystalDynamics.py#L385)

Questo metodo viene eseguito automaticamente SUBITO DOPO l'__init__.
Usato per calcolare le statistiche derivate (medie e dev st).

#### Signature

```python
def __post_init__(self): ...
```

### CrystalDynamicsResult().summary

[Show source in CrystalDynamics.py:403](../../zosojack/CMS/MolecularDynamics/CrystalDynamics.py#L403)

Stampa un riepilogo dei risultati della simulazione.

#### Signature

```python
def summary(self) -> None: ...
```