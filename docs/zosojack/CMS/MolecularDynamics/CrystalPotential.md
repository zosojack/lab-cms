# Crystalpotential

[CMS Index](../../../README.md#cms-index) / `zosojack` / [Cms](../index.md#cms) / [Moleculardynamics](./index.md#moleculardynamics) / Crystalpotential

> Auto-generated documentation for [zosojack.CMS.MolecularDynamics.CrystalPotential](../../../../zosojack/CMS/MolecularDynamics/CrystalPotential.py) module.

- [Crystalpotential](#crystalpotential)
  - [CrystalPotential](#crystalpotential)
    - [CrystalPotential().compute_forces](#crystalpotential()compute_forces)
    - [CrystalPotential().compute_potential](#crystalpotential()compute_potential)

## CrystalPotential

[Show source in CrystalPotential.py:122](../../../../zosojack/CMS/MolecularDynamics/CrystalPotential.py#L122)

Classe per calcolare il potenziale e le forze in una struttura cristallina. Implementa
il potenziale di Lennard-Jones per i primi vicini e un potenziale polinomiale di settimo
grado per i secondi vicini.

Attributes
----------
crystal : CrystalStructure
    Struttura cristallina di input.
sigma : float
    Parametro sigma del potenziale di Lennard-Jones (default: 2.644 Å).
epsilon : float
    Profondità del pozzo del potenziale di Lennard-Jones (default: 0.345 eV).
poly7 : PolynomialJunction | None
    Polinomio di giunzione di ordine 7 per le seconde vicine (opzionale).

Methods
-------
compute_potential() -> float
    Calcola il potenziale totale del cristallo.
compute_forces() -> np.ndarray
    Calcola le forze sugli atomi del cristallo.

#### Signature

```python
class CrystalPotential:
    def __init__(
        self,
        crystal: CrystalStructure,
        sigma: float = 2.644,
        epsilon: float = 0.345,
        poly7: PolynomialJunction = None,
    ) -> None: ...
```

### CrystalPotential().compute_forces

[Show source in CrystalPotential.py:188](../../../../zosojack/CMS/MolecularDynamics/CrystalPotential.py#L188)

Calcola le forze sugli atomi del cristallo.

#### Signature

```python
def compute_forces(self) -> np.ndarray: ...
```

### CrystalPotential().compute_potential

[Show source in CrystalPotential.py:158](../../../../zosojack/CMS/MolecularDynamics/CrystalPotential.py#L158)

Calcola il potenziale totale del cristallo.

#### Signature

```python
def compute_potential(self) -> float: ...
```