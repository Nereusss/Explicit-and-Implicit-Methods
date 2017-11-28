# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 15:47:20 2017

@author: Xingyun
"""

import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


def bs_finite_imp_call():
    
    print()
    print("Black-Scholes PDE, Numerical Estimation, Explicit Finite Differences")
    print()

    S_curr = 105.99
    K = 105.00
    sigma = 0.3128
    r = 0.0017
    T= 351.0/365.0
    divYield = 0.0170
    
    J = 200
    
    M = 15000
    
    S_max = 4*S_curr
    deltaS = S_max/J
    deltaT = T/M
    
    print("Strike Price K = %s" %K)
    print("Volatility sigma = %s" %sigma)
    print("Risk-free Rate r = %s" %r)
    print("Time to expiration T = %s" %T)
    print()
    
    print("Asset Price Steps J= %s" %J)
    print("Time Incriment Steps M = %s" %M)
    print("Size of Asset Price Steps = %s" %deltaS)
    print("Size of Time Steps = %s" %deltaT)
    
    print()
    
    B=np.zeros((J-1,J-1),dtype=float)
    S=np.zeros((J-1,1),dtype=float)
    C_hat=np.zeros((J-1,1),dtype=float)
   
    
    S[0]=deltaS
     
    j=1
    
    for j in range(1,J-1):
        S[j] = S[j-1]+deltaS
    
         
    j=0
    
    for j in range (0,J-1):
        C_hat[j] = max(S[j]-K,0)
        
           
    j=0
    k=0
    
    for j in range (0,J-1):
        for k in range (0,J-1):
            if k == j:
                B[j,j] = 1+(sigma**2)*((j+1)**2)*deltaT+r*deltaT
            elif k == j-1:
                B[j,k] = -0.5*((sigma**2)*((j+1)**2)*deltaT-(r-divYield)*(j+1)*deltaT)
            elif k == j+1:
                B[j,k] = -0.5*((sigma**2)*((j+1)**2)*deltaT+(r-divYield)*(j+1)*deltaT)
               
    Binverse = np.linalg.matrix_power(B,-1)    
           
    C_0 = 0
    C_start = C_hat
    
    
    m = 1
    
    while m <= M:
        
        C_1 = C_hat[0]
        C_Jminus1 = C_hat[J-2]
        
        C_hat = Binverse.dot(C_hat)
        
        C_hat[0] = (C_1 + 0.5*((sigma**2)*((1)**2)*deltaT-(r-divYield)*(1)*deltaT)*C_0+0.5*((sigma**2)*((1)**2)*deltaT+(r-divYield)*(1)*deltaT)*C_hat[1])/(1+(sigma**2)*((1)**2)*deltaT+r*deltaT)
        C_hat[J-2] = (C_Jminus1 + 0.5*((sigma**2)*((J-1)**2)*deltaT-(r-divYield)*(J-1)*deltaT)*C_hat[J-3]+0.5*((sigma**2)*((J-1)**2)*deltaT+(r-divYield)*(J-1)*deltaT)*(S_max-K*math.exp(-r*m*deltaT)))/(1+(sigma**2)*((J-1)**2)*deltaT+r*deltaT)
    
    
        if m == M/4:
           C1 = C_hat
        elif m == M/2:
           C2 = C_hat
        elif m == 3*M/4:
           C3 = C_hat
        elif m == M:
           C4 = C_hat
        
        
       
        m = m+1
    
    print("Asset Price = %s" %S[49])
    print("Call  Price = %s" %C_hat[49])

    
    plt.plot(S,C_start,'b--',S,C1,'r--',S,C2,'g--',S,C3,'r-',S,C4,'b--')
    
bs_finite_imp_call()