# PolynomialJunction.py
# Polynomial junction of order 7 between Lennard-Jones potential and zero
import numpy as np
class PolynomialJunction:
  """
  Calcola i coefficienti del polinomio di giunzione di ordine 7
  tra il potenziale di Lennard-Jones e lo zero, in funzione di
  R_P (punto di giunzione) e R_C (distanza di taglio).
  - continuità del potenziale in R_P e R_C
      • P7(R_P) = LJ(R_P)
      • P7(R_C) = 0
  - continuità della derivata prima [forze] in R_P e R_C 
  ...
  Restituisce i coefficienti A, B, C, D, E, F, G, H.
  """
        
  def __init__(self,
               R_C: float,
               R_P: float,
               epsilon: float = 0.345,
               sigma: float =  2.644):

    self.A = ( (1/((R_C - R_P)**7* R_P**12))*4* epsilon* R_C**4* sigma**6*
                 (2*R_P**6 *(-42* R_C**3 + 182*R_C**2* R_P -273* R_C* R_P**2 +
          143* R_P**3) + (455* R_C**3 - 1729*R_C**2* R_P +
          2223* R_C* R_P**2 - 969* R_P**3)* sigma**6) )
    
    self.B = ( (1/((R_C - R_P)**7* R_P**13))*16* epsilon* R_C**3* sigma**6* 
                (R_P**6* (54* R_C**4 - 154* R_C**3* R_P + 
      351* R_C *R_P**3 - 286* R_P**4) + (-315* R_C**4 + 749 *R_C**3 * R_P + 
      171 *R_C**2* R_P**2 - 1539* R_C* R_P**3 + 969* R_P**4)* sigma**6) )
                                            
    self.C = ( (1/((R_C - R_P)**7* R_P**14))*12* epsilon* R_C**2* sigma**6* 
      (R_P**6* (-63* R_C**5 - 7* R_C**4 *R_P + 665* R_C**3 *R_P**2 - 
      975* R_C**2* R_P**3 - 52* R_C* R_P**4 + 572* R_P**5) + 
      2 *(195* R_C**5 + 91* R_C**4* R_P - 1781* R_C**3* R_P**2 + 
      1995 *R_C**2* R_P**3 + 399* R_C* R_P**4 - 969* R_P**5)* sigma**6) )
          
    self.D = ( (1/((R_C - R_P)**7* R_P**15))*16* epsilon* sigma**6* 
      (R_C* R_P**6* (14* R_C**6 + 126* R_C**5* R_P - 420* R_C**4* R_P**2 
      -90* R_C**3* R_P**3 + 1105* R_C**2* R_P**4 - 624* R_C* R_P**5 - 
      286 *R_P**6) + R_C* (-91* R_C**6 - 819* R_C**5* R_P + 2145 * 
      R_C**4 * R_P**2 + 1125* R_C**3* R_P**3 - 5035* R_C**2* R_P**4 + 
      1881* R_C* R_P**5 + 969* R_P**6)* sigma**6) )
                              
    self.E = ( (1/((R_C - R_P)**7* R_P**15))*4* epsilon* sigma**6* 
      (2* R_P**6* (-112* R_C**6 - 63* R_C**5* R_P + 1305* R_C**4* R_P**2 
      -1625* R_C**3* R_P**3 - 585* R_C**2* R_P**4 + 
      1287 *R_C* R_P**5 + 143* R_P**6) + (1456*R_C**6 +1404*R_C**5* R_P - 
      14580 *R_C**4* R_P**2 + 13015* R_C**3* R_P**3 + 7695* R_C**2* R_P**4 - 
      8721 *R_C* R_P**5 - 969* R_P**6)* sigma**6) )
                                              
    self.F = ( (1/((R_C - R_P)**7* R_P**15))*48* epsilon* sigma**6* 
      (-R_P**6* (-28* R_C**5 + 63* R_C**4* R_P + 65* R_C**3* R_P**2 - 
      247* R_C**2* R_P**3 + 117* R_C* R_P**4 + 65* R_P**5) + 
      (-182* R_C**5 + 312* R_C**4* R_P + 475* R_C**3* R_P**2 - 
      1140* R_C**2* R_P**3 + 342* R_C* R_P**4 + 228* R_P**5)* sigma**6) )
          
    self.G = ( (1/((R_C - R_P)**7* R_P**15))*4* epsilon* sigma**6* (R_P**6* 
      (-224* R_C**4 + 819* R_C**3* R_P - 741* R_C**2* R_P**2 
      -429* R_C* R_P**3 + 715* R_P**4) + 2 *(728* R_C**4 
      -2223* R_C**3* R_P + 1425* R_C**2* R_P**2 
      +1292* R_C* R_P**3 - 1292* R_P**4)* sigma**6) )
                          
    self.H = ( (1/((R_C - R_P)**7* R_P**15))*16* epsilon* sigma**6* (R_P**6* 
      (14* R_C**3 - 63* R_C**2* R_P + 99* R_C* R_P**2 - 55* R_P**3) + 
      (-91* R_C**3 + 351* R_C**2* R_P - 459* R_C* R_P**2 
      +204* R_P**3)* sigma**6) )
  
  @property
  def coeffs_array(self):
    return np.array([self.A, self.B, self.C, self.D, self.E, self.F, self.G, self.H])
  
  def eval(self, r):
    return (self.A + self.B*r + self.C*r**2 + self.D*r**3 + 
    self.E*r**4 + self.F*r**5 + self.G*r**6 + self.H*r**7)
    
  def eval_derivative(self, r):
    return (self.B + 2*self.C*r + 3*self.D*r**2 + 4*self.E*r**3 + 
    5*self.F*r**4 + 6*self.G*r**5 + 7*self.H*r**6)