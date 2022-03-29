# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 16:02:50 2022

@author: ellen
"""
#testing stability of steady state
#is y_eqm steady?

import numpy as np
import sympy as sym
from math import comb
n=2
m=3

##y_0 that gives all positive steady states (initial y)
#y=np.array([[1.91000000e+02],
#       [6.22000000e+02],
#       [5.24000000e+02],
#       [8.87000000e+02],
#       [9.23000000e+02],
#       [6.65202821e+05],
#       [4.79546424e+05],
#       [7.97864685e+05]])

##those resulting steady states
#y_eqm=np.array([[0.29167146],
#       [0.05132694],
#       [0.88146805],
#       [1.02758093],
#       [0.04774258],
#       [4.82453013],
#       [4.48961497],
#       [4.93698289]])


T_1 = sym.Symbol('T_1')
T_2 = sym.Symbol('T_2')
P_1=sym.Symbol('P_1')
P_2=sym.Symbol('P_2')
P_3=sym.Symbol('P_3')
C_12=sym.Symbol('C_12')
C_13=sym.Symbol('C_13')
C_23=sym.Symbol('C_23')

y_s=np.array([[T_1],[T_2],[P_1],[P_2],[P_3],[C_12],[C_13],[C_23]])

#loops wont allow the symbol entries
#r(y_s, sigma, A, f,)
#check the eigenvalues of the A part, but where does the specfic y_eqm values come into this?
#if we take the eigenvalues of when it also multiplies by the y_eqm?
A_eig=np.linalg.eigvals(A)

#maybe could eventually initialise system such that the real parts of these eigenvalues are all neg
if ODE==0:
    ## Jacobian of g
    G=np.zeros([comb(m,2),m])
    for i in range(0,m):
        for k in range(i+1,m):
            if i==0:
                l=(i)*(m-i)+k-i-1
            else:
                l=l+1
            print(l)
            G[l,i]=rho_vec[l]*y_eqm[n+k]
            G[l,k]=rho_vec[l]*y_eqm[n+i]


if ODE==1:
    ## Jacobian of g
    #G=np.zeros([comb(m,2),m], dtype=object) #so it will take symbols
    G=np.zeros([comb(m,2),m])
    for i in range(0,m):
        for k in range(i+1,m):
            if i==0:
                l=(i)*(m-i)+k-i-1
            else:
                l=l+1
            #print(l)
            G[l,i]=y_eqm[n+1] # or use y_s
            G[l,k]=y_eqm[n+k]
        
#entire Jacobian of g        
#D_g=np.zeros([n+m+comb(m,2),n+m+comb(m,2)],dtype=object)
D_g=np.zeros([n+m+comb(m,2),n+m+comb(m,2)])
D_g[n+m:,n:n+m]=G

#eigenvalues of entire Jacobian of system #f+Ay+g(y) assesed at y_eqm
eqm_eigs=np.linalg.eigvals(A+D_g)
print(eqm_eigs)






