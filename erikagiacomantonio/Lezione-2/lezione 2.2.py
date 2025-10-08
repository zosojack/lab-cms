#plot ept /N vs N
import numpy as np
import matplotlib.pyplot as plt

file = "Epot.txt"
graph = np.loadtxt(file)

x = graph[:, 0]
y = graph[:, 2]

plt.plot(x, y, marker="o", label= "senza cutoff")


file2 = "Epot_cutoff.txt"
graph_cutoff = np.loadtxt(file2)
x2 = graph_cutoff[:, 0]
y2 = graph_cutoff[:, 2]
plt.plot(x2, y2, marker="o", label = "con cutoff")
plt.xlabel("Natom")
plt.ylabel("Epot/N")
plt.legend()
plt.title("Energia potenziale per atomo in funzione del numero di atomi, con r_c= 3")
plt.show()

# Epot ha questo andamento perché all'aumentare del numero di atomi ho più legami, col cutoff si alzerà 
