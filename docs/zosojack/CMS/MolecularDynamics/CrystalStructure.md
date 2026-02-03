# Crystalstructure

[CMS Index](../../../README.md#cms-index) / `zosojack` / [Cms](../index.md#cms) / [Moleculardynamics](./index.md#moleculardynamics) / Crystalstructure

> Auto-generated documentation for [zosojack.CMS.MolecularDynamics.CrystalStructure](../../../../zosojack/CMS/MolecularDynamics/CrystalStructure.py) module.

- [Crystalstructure](#crystalstructure)
  - [CrystalStructure](#crystalstructure)
    - [CrystalStructure().add_atom](#crystalstructure()add_atom)
    - [CrystalStructure().copy](#crystalstructure()copy)
    - [CrystalStructure().crystal_center](#crystalstructure()crystal_center)
    - [CrystalStructure().displacements_matrix](#crystalstructure()displacements_matrix)
    - [CrystalStructure.empty](#crystalstructureempty)
    - [CrystalStructure().find_neighbours](#crystalstructure()find_neighbours)
    - [CrystalStructure.from_file](#crystalstructurefrom_file)
    - [CrystalStructure().only_neighbours_distance](#crystalstructure()only_neighbours_distance)
    - [CrystalStructure().print_neighbours](#crystalstructure()print_neighbours)
    - [CrystalStructure().print_second_neighbours](#crystalstructure()print_second_neighbours)
    - [CrystalStructure().set_R_C](#crystalstructure()set_r_c)
    - [CrystalStructure().set_R_P](#crystalstructure()set_r_p)
    - [CrystalStructure().set_R_V](#crystalstructure()set_r_v)
    - [CrystalStructure().set_pbc](#crystalstructure()set_pbc)
    - [CrystalStructure().vec_x](#crystalstructure()vec_x)
    - [CrystalStructure().vec_y](#crystalstructure()vec_y)
    - [CrystalStructure().vec_z](#crystalstructure()vec_z)
  - [_find_neighbour_masks_kernel](#_find_neighbour_masks_kernel)
  - [_only_neighbours_distance_kernel](#_only_neighbours_distance_kernel)

## CrystalStructure

[Show source in CrystalStructure.py:180](../../../../zosojack/CMS/MolecularDynamics/CrystalStructure.py#L180)

CrystalStructure
================
Classe per rappresentare una struttura cristallina.

Attributes
----------
positions : np.ndarray
    Matrice Nx3 delle posizioni atomiche.
N_atoms : int
    Numero totale di atomi.
R_C : float
    Distanza di taglio per i primi vicini.
R_P : float
    Punto di giunzione polinomiale che separa primi e secondi vicini.
R_V : float
    Raggio della Verlet cage entro cui tenere traccia degli atomi.
neighbours_computed : bool
    Flag che indica se matrici di vicini e distanze sono aggiornate.
reference_positions : np.ndarray | None
    Copia delle posizioni quando i vicini sono stati calcolati.
pbc : np.ndarray
    Dimensioni della cella per le condizioni periodiche.

Properties
----------
vec_x : np.ndarray
    Vettore delle coordinate x degli atomi.
vec_y : np.ndarray
    Vettore delle coordinate y degli atomi.
vec_z : np.ndarray
    Vettore delle coordinate z degli atomi.
displacements_matrix : np.ndarray
    Tensore NxNx3 dei vettori spostamento fra tutti gli atomi.
crystal_center : np.ndarray
    Coordinate (x, y, z) del centro del volume del cristallo.

Methods
-------
from_file(filename) -> CrystalStructure
    Crea un oggetto CrystalStructure leggendo da file.
empty(n_atoms) -> CrystalStructure
    Costruttore alternativo: crea un oggetto CrystalStructure vuoto.
copy() -> CrystalStructure
    Restituisce una copia di un oggetto CrystalStructure.
set_R_C(R_C) -> None
    Imposta la distanza di taglio R_C.
set_R_P(R_P) -> None
    Imposta il punto di giunzione polinomiale R_P.
set_R_V(R_V) -> None
    Imposta la grandezza della Verlet cage R_V.
set_pbc(pbc) -> None
    Imposta la dimensione della cella per la periodicità al contorno.
add_atom(position) -> None
    Aggiunge un atomo alla struttura cristallina.
find_neighbours() -> None
    Ricalcola tutte le distanze e trova primi e secondi vicini.
only_neighbours_distance() -> None
    Ricalcola solo le distanze dai vicini già noti.
print_neighbours(index=None) -> None
    Stampa gli indici dei primi vicini per ogni atomo o di uno specifico.
print_second_neighbours(index=None) -> None
    Stampa gli indici dei secondi vicini per ogni atomo o di uno specifico.

#### Signature

```python
class CrystalStructure:
    def __init__(self, positions): ...
```

### CrystalStructure().add_atom

[Show source in CrystalStructure.py:420](../../../../zosojack/CMS/MolecularDynamics/CrystalStructure.py#L420)

Aggiunge un atomo alla struttura cristallina.

Parameters
----------
position : array-like
    Coordinate (x, y, z) dell'atomo da aggiungere.

#### Signature

```python
def add_atom(self, position) -> None: ...
```

### CrystalStructure().copy

[Show source in CrystalStructure.py:354](../../../../zosojack/CMS/MolecularDynamics/CrystalStructure.py#L354)

#### Signature

```python
def copy(self) -> CrystalStructure: ...
```

### CrystalStructure().crystal_center

[Show source in CrystalStructure.py:326](../../../../zosojack/CMS/MolecularDynamics/CrystalStructure.py#L326)

Restituisce le coordinate (x, y, z) del centro del volume del cristallo.

#### Signature

```python
@property
def crystal_center(self) -> np.ndarray: ...
```

### CrystalStructure().displacements_matrix

[Show source in CrystalStructure.py:295](../../../../zosojack/CMS/MolecularDynamics/CrystalStructure.py#L295)

Restituisce il tensore NxNx3 dei vettori spostamento fra tutti gli atomi.
Shape: (N_atomi, N_atomi, 3)

#### Signature

```python
@property
def displacements_matrix(self) -> np.ndarray: ...
```

### CrystalStructure.empty

[Show source in CrystalStructure.py:270](../../../../zosojack/CMS/MolecularDynamics/CrystalStructure.py#L270)

Costruttore alternativo: crea Crystal vuoto

#### Signature

```python
@classmethod
def empty(cls, n_atoms) -> CrystalStructure: ...
```

### CrystalStructure().find_neighbours

[Show source in CrystalStructure.py:436](../../../../zosojack/CMS/MolecularDynamics/CrystalStructure.py#L436)

Ricalcola OGNI distanza e trova primi e secondi vicini per ogni atomo in base a R_C e R_P.
Gli atomi entro una distanza di Verlet (R_C < r < R_V) sono aggiunti alla matrice dei secondi vicini,
ma la loro distanza non è salvata in distance_matrix.

Returns
-------
None

#### Signature

```python
def find_neighbours(self) -> None: ...
```

### CrystalStructure.from_file

[Show source in CrystalStructure.py:262](../../../../zosojack/CMS/MolecularDynamics/CrystalStructure.py#L262)

Crea un Crystal leggendo da file.

#### Signature

```python
@classmethod
def from_file(cls, filename) -> CrystalStructure: ...
```

### CrystalStructure().only_neighbours_distance

[Show source in CrystalStructure.py:460](../../../../zosojack/CMS/MolecularDynamics/CrystalStructure.py#L460)

Ricalcola solamente le distanze dai vicini già noti, senza aggiornarli.
Se le posizioni sono cambiate troppo, neighbours_computed è posto False ed
è necessario rieseguire find_neighbours().

Returns
-------
None

#### Signature

```python
def only_neighbours_distance(self) -> None: ...
```

### CrystalStructure().print_neighbours

[Show source in CrystalStructure.py:483](../../../../zosojack/CMS/MolecularDynamics/CrystalStructure.py#L483)

Stampa gli indici dei primi vicini per ogni atomo o di uno nello specifico.

#### Signature

```python
def print_neighbours(self, index=None) -> None: ...
```

### CrystalStructure().print_second_neighbours

[Show source in CrystalStructure.py:501](../../../../zosojack/CMS/MolecularDynamics/CrystalStructure.py#L501)

Stampa gli indici dei secondi vicini per ogni atomo o di uno nello specifico.

#### Signature

```python
def print_second_neighbours(self, index=None) -> None: ...
```

### CrystalStructure().set_R_C

[Show source in CrystalStructure.py:357](../../../../zosojack/CMS/MolecularDynamics/CrystalStructure.py#L357)

Imposta la distanza di taglio R_C usata per trovare i vicini.

Parameters
----------
R_C : float
    Distanza di taglio per i primi vicini.

#### Signature

```python
def set_R_C(self, R_C) -> None: ...
```

### CrystalStructure().set_R_P

[Show source in CrystalStructure.py:368](../../../../zosojack/CMS/MolecularDynamics/CrystalStructure.py#L368)

Imposta il punto di giunzione polinomiale R_P;
divide primi e 'secondi' vicini.

Parameters
----------
R_P : float
    Raggio che separa primi e secondi vicini.

#### Signature

```python
def set_R_P(self, R_P) -> None: ...
```

### CrystalStructure().set_R_V

[Show source in CrystalStructure.py:380](../../../../zosojack/CMS/MolecularDynamics/CrystalStructure.py#L380)

Imposta la grandezza della Verlet cage R_V.
Consigliato: R_V = R_C + 0.5 Å.
Non può essere minore di R_C, altrimenti solleva un ValueError.

Parameters
----------
R_V : float
    Raggio della Verlet cage (deve essere >= R_C).

#### Signature

```python
def set_R_V(self, R_V) -> None: ...
```

### CrystalStructure().set_pbc

[Show source in CrystalStructure.py:395](../../../../zosojack/CMS/MolecularDynamics/CrystalStructure.py#L395)

Imposta la dimensione della cella per la condizione di periodicità al contorno.
Deve essere un array-like di 3 elementi.
Deve essere consistente con R_C.
Può ricevere np.inf per indicare nessuna periodicità in una direzione.

Parameters
----------
pbc : array-like
    Dimensioni della cella (Lx, Ly, Lz) in Å; usare np.inf per direzione non periodica.

#### Signature

```python
def set_pbc(self, pbc) -> None: ...
```

### CrystalStructure().vec_x

[Show source in CrystalStructure.py:280](../../../../zosojack/CMS/MolecularDynamics/CrystalStructure.py#L280)

Restituisce il vettore delle coordinate x degli atomi.

#### Signature

```python
@property
def vec_x(self) -> np.ndarray: ...
```

### CrystalStructure().vec_y

[Show source in CrystalStructure.py:285](../../../../zosojack/CMS/MolecularDynamics/CrystalStructure.py#L285)

Restituisce il vettore delle coordinate y degli atomi.

#### Signature

```python
@property
def vec_y(self) -> np.ndarray: ...
```

### CrystalStructure().vec_z

[Show source in CrystalStructure.py:290](../../../../zosojack/CMS/MolecularDynamics/CrystalStructure.py#L290)

Restituisce il vettore delle coordinate z degli atomi.

#### Signature

```python
@property
def vec_z(self) -> np.ndarray: ...
```



## _find_neighbour_masks_kernel

[Show source in CrystalStructure.py:47](../../../../zosojack/CMS/MolecularDynamics/CrystalStructure.py#L47)

HACK: poiché rij viene salvato solo per i vicini, i confronti vengono effettuati
senza calcolarne la radice quadrata, per efficienza. Per farlo, anche R_P, R_C e R_V
vanno considerati al quadrato.

#### Signature

```python
@njit(cache=True)
def _find_neighbour_masks_kernel(
    positions, R_P, R_C, R_V, pbc
) -> tuple[np.ndarray, np.ndarray, np.ndarray]: ...
```



## _only_neighbours_distance_kernel

[Show source in CrystalStructure.py:105](../../../../zosojack/CMS/MolecularDynamics/CrystalStructure.py#L105)

#### Signature

```python
@njit(cache=True)
def _only_neighbours_distance_kernel(
    positions, R_P, R_C, R_V, neighbour, second, pbc
) -> tuple[np.ndarray, np.ndarray, np.ndarray]: ...
```