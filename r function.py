# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 13:41:30 2022

@author: ellen
"""

import numpy as np
from math import comb
from scipy import optimize

#can I put my whole g calculation inside of this function?
n=2
m=3
ODE=1
#rho=np.zeros([m,m])


#def r(y):
#    for i in range(m):
#        for j in range(m):
#            if j >i:
#                rho[i,j]=sigma[i,j]*(y[i+n+1]/y[j+n+1])+sigma[j,i]*(y[i+n+1]/y[i+n+1])
#return np.array(f+np.dot(A,x)+g)

#function to calculate the g vector per every y
def r(y, sigma, A, f):
    if ODE==0:
        rho=np.zeros([m,m])
        for i in range(m): #upper bound value not included in range of loop
            for j in range(m):
                if j >i:
                    rho[i,j]=sigma[i,j]*(y[i+n]/y[j+n])+sigma[j,i]*(y[j+n]/y[i+n])
                    #dont forget about 0'th element 
        rho_vec=[]        
        for i in range(m):
            for j in range(m):
                if rho[i,j] !=0:
                    rho_entry=rho[i,j]
                    rho_vec.append(rho_entry)
        rho_vec=np.transpose(np.asarray([rho_vec])) 
        
        PP_list=[]
        for i in range(m):
            for j in range(m):
                if rho[i,j] !=0: 
                    PP_entry=y[i+n]*y[j+n]
                    PP_list.append(PP_entry)          
        PP=np.zeros([m,m])
        for i in range(m):
            for j in range(m):
                PP[i,i]=PP_list[i] 
        g_t=np.zeros([n+m,1])
        g_b=np.dot(PP,rho_vec)
        g=np.vstack((g_t,g_b))
    
    
    if ODE==1:
        rho=np.ones([m,m])
        rho_vec=np.ones([comb(m,2),1])
    
    PP=[]
    for i in range(m):
        for j in range(i+1,m):
            PP.append(y[i+n]*y[j+n])
    PP=np.array(PP) 
    
    g_t=np.zeros([n+m,1])
    g_b=PP
    g=np.vstack((g_t,g_b))
    return np.array(f+np.dot(A,y)+g)
#so we have rho[0,1]*P[0]*P[1] as rho_12 P_1 P_2 due to 0'th element

#solves system of equations defined in function for being equal to zero, stea
#arguments defined in Basic pathogen combination model file
y_eqm = optimize.newton(r, 10*np.ones([8,1]),args=(sigma, A, f),tol=1.48e-2)
#works even with default tolerance                
