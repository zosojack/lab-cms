n=10
for i in range(1,n+1):
    print("value of i:",i)
    print("value")

print("2value of i:", i) #l'indentazione è impo

print("------------------------")

mysum = 0 #definisci sempre prima le variabili
for i in range(1,n+1):
    print("value of i:",i)
    mysum = mysum + i
print("Sum value is:",mysum)

print("------------------------")

#with use of a built-in function (sum) of Phython standard without using external libraries
mysum2 = 0
mysum2 = sum(range(1,n+1))
print("sum value is:",mysum2)

print("------------------------")

mysum3 = 0
for i in range(1, n+1):
    mysum3 += i #scritto così te lo somma ogni volta di 1
    if mysum3 > 7:
        print("partial value of the sum",mysum3)

print("------------------------")
mysum4 = 0
for i in range(1, n+1):
    mysum4 += i 
    if mysum4 > 3 and mysum4 < 14:
        print("now",i, mysum4)
        
print("------------------------")
mysum5 = 0
for i in range(1, n+1):
    mysum5 += i 
    if mysum5 > 7:
        print("partial value of the sum:",mysum5)
        break
print("final value of the sum:", mysum5)

print("------------------------")
mysum6 = 0
for i in range(1, n+1):
    mysum6 += i 
    if mysum6 > 7:
        print("greater than 7")
    else:
        print("not greater")
        
print("final value of the sum:", mysum6)
    
print("------------------------")   
mysum7 = 0
for i in range(1, n+1):
    mysum7 += i 
    if mysum7 > 7:
        print("greater than 7")
    elif mysum7 > 2:
        print("between 7 and 2") #•alternativa all'if ma non else
    else:
        print("not greater")
        
print("final value of the sum:", mysum7)   
    
print("------------------------")
import numpy as np   
filename = "fcc100a108.txt"
#con spider puoi aprire i file di testo in questa finestra
pos = np.loadtxt(filename)
print("total numer of atoms", pos.shape[0]) #indico solo il numero di righe e non di colonne
print("Y coordinate fo atmo 54:", pos[53,1]) #riga, colonna
natoms= pos.shape[0] 
#dx= difference in the x coordinate of atom 3 and of atom 4
dx = pos[2,0]-pos[3,0]
print("diff pos x between atom 3 and 4:", dx)

#label of atom i label of atom j distance between atom i and j
for i in range(1, natoms):
    for j in range(1, natoms):
        dx = pos[i,0]-pos[j,0]
        dy = pos[i,1]-pos[j,1]
        dz = pos[i,2]-pos[j,2]
        d2 = dx*dx + dy*dy + dz*dz
        d2 = d2**0.5
        print("atom", i, " atom", j," distance is:" , d2)
        
print("------------------------")
        
#trova la distanza minima, e ci serve perché sarà il primo vicino
#ovviamente inserisci il caso in cui la distanza non può essere 0 perché non vale considerare se stessi i \neq 0

#prima utilizzo una prima guess, metto un numero molto grande
dmin = 6326551545   
for i in range(1, natoms):
    for j in range(1, natoms):
        dx = pos[i,0]-pos[j,0]
        dy = pos[i,1]-pos[j,1]
        dz = pos[i,2]-pos[j,2]
        d2 = dx*dx + dy*dy + dz*dz
        d2 = d2**0.5
        if d2 < dmin and i !=j:
            dmin = d2
        
print("dmin" , dmin) #questa distanza sarà il mio 1st nn

#trova anche il 2nn e il 3 nn e la distanza massima 
    
    
    
    
    
    
    
    
    