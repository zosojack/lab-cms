filename_list = ["fcc100a108.txt", "fcc100a256.txt", "fcc100a500.txt", "fcc100a864.txt", "fcc100a1372.txt", "fcc100a2048.txt"]

import numpy as np
for filename_list in filename_list:
    pos = np.loadtxt(filename_list)
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
            if i!=j and d2 < 3:
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

