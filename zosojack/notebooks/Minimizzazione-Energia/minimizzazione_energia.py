from pathlib import Path
import sys

# aggiungi la root (lab-cms) a sys.path
ROOT = Path(__file__).resolve().parents[3]  # /Users/zosojack/lab-cms
sys.path.insert(0, str(ROOT))

from zosojack.libraries.CrystalStructure import CrystalStructure as Crystal
from zosojack.libraries.SteepestDescend import SteepestDescend

# Nome del file (numero di atomi da studiare: 256)
filename = 'zosojack/data/fcc100a256.txt'

# interrompe l'algoritmo se non converge entro max_steps
max_steps = 5000

# coefficiente di discesa
C_steep = 0.005

# forza limite per la convergenza
F_tol = 1e-3 

# Inizializza la struttura cristallina
cristallo = Crystal.from_file(filename)
steepest = SteepestDescend(cristallo)
steepest.minimize_energy(C_steep=C_steep, F_tol=F_tol, max_steps=max_steps, pol_junction=False)
