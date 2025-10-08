# leggi le posizioni atomiche dal file fcc100a108
#usa E pot = 1/2 \sum:i,j=1...N phi(r_ij) con phi lennard jones potential 
# calcoliamo l'energia potenziale totale di un FCC 
# rifallo con tutte le celle disponibili 256 500 864 1372 2048
# fai un file di testo con l'editor di spider dove metti 108 epot ed epot/108 CHE è L'EPOT per atomo così per tutti i numeri 
# poi plotta l'energia potenziale per atomo vs N
#ripeti tutto con il cut of rc = 3
#stampa per quella da 108 atom i, forza lungo x(i) y(i) z(i) sempre col cutoff

import numpy as np
filename = "fcc100a108.txt"
pos = np.loadtxt(filename)
print("total numer of atoms", pos.shape[0]) #indico solo il numero di righe e non di colonne
natoms= pos.shape[0] 

#calcolo le distanze tra gli atomi per i diverso da j e creo l'array contenente le distanze
d2vec = []
for i in range(natoms):
    for j in range(natoms):
        dx = pos[i,0]-pos[j,0]
        dy = pos[i,1]-pos[j,1]
        dz = pos[i,2]-pos[j,2]
        d2 = (dx*dx) + (dy*dy) + (dz*dz)
        d2 = d2**0.5
        if i!=j:
            d2vec.append(d2)
            
d2vec = np.array(d2vec)
ndistances = d2vec.shape[0]

#ora calcolo il potenziale di LJ per ogni posizione
eps = 0.345 #ev
sig = 2.644 #A
phirij = 0
phirijvec = []
for i in range(ndistances):
    phirij = 4*eps*((sig/d2vec[i])**12 - (sig/d2vec[i])**6)
    phirijvec.append(phirij)
    
phirijvec = np.array(phirijvec)
nphi = phirijvec.shape[0]

#sommo per calcolare il potenziale totale
Epot_tot = 0.5*np.sum(phirijvec)
Epot_peratomo = Epot_tot/natoms
print("epot tot", Epot_tot)
print("epot per atomo", Epot_peratomo)

print("----------------------")
filename2 = "fcc100a256.txt"
pos = np.loadtxt(filename2)
print("total numer of atoms", pos.shape[0]) #indico solo il numero di righe e non di colonne
natoms2= pos.shape[0] 