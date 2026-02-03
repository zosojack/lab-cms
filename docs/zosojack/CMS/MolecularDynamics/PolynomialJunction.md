# Polynomialjunction

[CMS Index](../../../README.md#cms-index) / `zosojack` / [Cms](../index.md#cms) / [Moleculardynamics](./index.md#moleculardynamics) / Polynomialjunction

> Auto-generated documentation for [zosojack.CMS.MolecularDynamics.PolynomialJunction](../../../../zosojack/CMS/MolecularDynamics/PolynomialJunction.py) module.

- [Polynomialjunction](#polynomialjunction)
  - [PolynomialJunction](#polynomialjunction)
    - [PolynomialJunction().coeffs_array](#polynomialjunction()coeffs_array)
    - [PolynomialJunction().eval](#polynomialjunction()eval)
    - [PolynomialJunction().eval_derivative](#polynomialjunction()eval_derivative)

## PolynomialJunction

[Show source in PolynomialJunction.py:4](../../../../zosojack/CMS/MolecularDynamics/PolynomialJunction.py#L4)

PolynomialJunction
=====================
Classe che implementa la giunzione polinomiale di ordine 7 tra il potenziale di Lennard-Jones e zero.
Calcola i coefficienti del polinomio di giunzione di ordine 7 tra il potenziale di Lennard-Jones e lo zero,
in funzione di R_P (punto di giunzione) e R_C (distanza di taglio).

Vincoli imposti:
----------------
- continuità del potenziale in R_P e R_C
    • P7(R_P) = LJ(R_P)
    • P7(R_C) = 0
- continuità della derivata prima [forze] in R_P e R_C
...
Restituisce i coefficienti A, B, C, D, E, F, G, H.

#### Signature

```python
class PolynomialJunction:
    def __init__(
        self, R_C: float, R_P: float, epsilon: float = 0.345, sigma: float = 2.644
    ): ...
```

### PolynomialJunction().coeffs_array

[Show source in PolynomialJunction.py:75](../../../../zosojack/CMS/MolecularDynamics/PolynomialJunction.py#L75)

#### Signature

```python
@property
def coeffs_array(self) -> np.ndarray: ...
```

### PolynomialJunction().eval

[Show source in PolynomialJunction.py:79](../../../../zosojack/CMS/MolecularDynamics/PolynomialJunction.py#L79)

#### Signature

```python
def eval(self, r) -> float: ...
```

### PolynomialJunction().eval_derivative

[Show source in PolynomialJunction.py:83](../../../../zosojack/CMS/MolecularDynamics/PolynomialJunction.py#L83)

#### Signature

```python
def eval_derivative(self, r) -> float: ...
```