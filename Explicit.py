# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 13:19:18 2017

@author: Xingyun
"""

import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


def bs_finite_exp_call_put():
    
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
    
    M = 4000
    
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
    
    A=np.zeros((J-1,J-1),dtype=float)
    S=np.zeros((J-1,1),dtype=float)
    C_hat_call=np.zeros((J-1,1),dtype=float)
    C_hat_put=np.zeros((J-1,1),dtype=float)
    
    S[0]=deltaS
     
    j=1
    
    for j in range(1,J-1):
        S[j] = S[j-1]+deltaS
    
         
    j=0
    
    for j in range (0,J-1):
        C_hat_call[j] = max(S[j]-K,0)
        C_hat_put[j] = max(K-S[j],0)
           
    j=0
    k=0
    
    for j in range (0,J-1):
        for k in range (0,J-1):
            if k == j:
                A[j,j] = 1-(sigma**2)*((j+1)**2)*deltaT-r*deltaT
            elif k == j-1:
                A[j,k] = 0.5*((sigma**2)*((j+1)**2)*deltaT-(r-divYield)*(j+1)*deltaT)
            elif k == j+1:
                A[j,k] = 0.5*((sigma**2)*((j+1)**2)*deltaT+(r-divYield)*(j+1)*deltaT)
               
                
    C_0 = 0
    C_J = 0
    C_start_call = C_hat_call
    C_start_put = C_hat_put
    
    m = 1
    
    while m <= M:
        C_min_call = 0.5*((sigma**2)*((1)**2)*deltaT-(r-divYield)*(1)*deltaT)*C_0+(1-(sigma**2)*((1)**2)*deltaT-r*deltaT)*C_hat_call[0]+0.5*((sigma**2)*((1)**2)*deltaT+(r-divYield)*(1)*deltaT)*C_hat_call[1]
        C_max_call = 0.5*((sigma**2)*((J-1)**2)*deltaT-(r-divYield)*(J-1)*deltaT)*C_hat_call[J-3]+(1-(sigma**2)*((J-1)**2)*deltaT-r*deltaT)*C_hat_call[J-2]+0.5*((sigma**2)*((J-1)**2)*deltaT+(r-divYield)*(J-1)*deltaT)*(S_max-K*math.exp(-r*m*deltaT))
        
        C_min_put = 0.5*((sigma**2)*((1)**2)*deltaT-(r-divYield)*(1)*deltaT)*(K*math.exp(-r*m*deltaT))+(1-(sigma**2)*((1)**2)*deltaT-r*deltaT)*C_hat_put[0]+0.5*((sigma**2)*((1)**2)*deltaT+(r-divYield)*(1)*deltaT)*C_hat_put[1]
        C_max_put = 0.5*((sigma**2)*((J-1)**2)*deltaT-(r-divYield)*(J-1)*deltaT)*C_hat_put[J-3]+(1-(sigma**2)*((J-1)**2)*deltaT-r*deltaT)*C_hat_put[J-2]+0.5*((sigma**2)*((J-1)**2)*deltaT+(r-divYield)*(J-1)*deltaT)*C_J
        
        C_hat_call = A.dot(C_hat_call)
        C_hat_put = A.dot(C_hat_put)
    
        C_hat_call[0] = C_min_call
        C_hat_call[J-2] = C_max_call
        C_hat_put[0] = C_min_put
        C_hat_put[J-2] = C_max_put          
                  
    
        if m == M/4:
           C1 = C_hat_call
        elif m == M/2:
           C2 = C_hat_call
        elif m == 3*M/4:
           C3 = C_hat_call
        elif m == M:
           C4 = C_hat_call
           
        if m == M/4:
           C5 = C_hat_put
        elif m == M/2:
           C6 = C_hat_put
        elif m == 3*M/4:
           C7 = C_hat_put
        elif m == M:
           C8 = C_hat_put
    
        m = m+1
    
    print("Asset Price = %s" %S[49])
    print("Call price = %s" %C_hat_call[49])
    print("Put price = %s" %C_hat_put[49])
    
    plt.plot(S,C_start_call,'b--',S,C1,'r--',S,C2,'g--',S,C3,'r-',S,C4,'b--',S,C_start_put,'b--',S,C5,'r--',S,C6,'g--',S,C7,'r-',S,C8,'b--')
    
bs_finite_exp_call_put()