# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 16:50:23 2022

@author: ellen
"""
import numpy as np
import random 
from math import comb
from scipy import optimize

n=2 #number of tree species
m=3 #number of pathogen species
T=np.transpose(np.array([random.sample(range(0, 1000), n)])) #amounts of each tree species
P=np.transpose(np.array([random.sample(range(400, 1000), m)])) #amount of each pathogen species
#what should the relationship be between the numbers of pathogen and tree speices?
#what is the order of magnitude difference?

## there needs to be less in each combo class than in each of the contributing
##pathogen classes
#current method of ensuring this is crude, need to fix, relate to the Ps
C=np.transpose(np.array([random.sample(range(0, 400), comb(m,2))]))
#create vector of all types
y= np.concatenate((T, P, C))


##death and birth rate vectors
a=np.random.random((n,1))
b=np.random.random((m,1))
d=np.random.random((comb(m,2),1))
f=np.concatenate((a,-b,-d))

##interaction constants and A matrix
alpha=np.random.random((n,m))
gamma=np.random.random((n,comb(m,2)))
beta=np.random.random((m,n))
sigma=np.random.random((m,comb(m,2)))
sigma[np.diag_indices_from(sigma)] = 0      #sets diagonal elements as zero
phi=np.random.random((comb(m,2),n))
A_r1=np.hstack((np.zeros([n,n]),-alpha,-gamma))
A_r2=np.hstack((beta,np.zeros([m,m]),-sigma))
A_r3=np.hstack((phi,np.zeros([comb(m,2),m]),np.zeros([comb(m,2),comb(m,2)])))           
A=np.vstack((A_r1,A_r2,A_r3))                

##create rho vector/matrix
#this incorporates a constraint to do with energy conservation
#between energy given by each pathogen to its interaction with another
#pathogen as they attack the tree
rho=np.zeros([m,m])
for i in range(m):
    for j in range(m):
        if j >i: #to count non-zero elements exactly once
            rho[i,j]=sigma[i,j]*(P[i]/P[j])+sigma[j,i]*(P[j]/P[i])
            
        
rho_vec=[]        
#rho_vec=np.zeros([m,1])
for i in range(m):
    for j in range(m):
        if rho[i,j] !=0:
            rho_entry=rho[i,j]
            rho_vec.append(rho_entry)
rho_vec=np.transpose(np.asarray([rho_vec]))

##create quadratic term matrix
#PP=np.zeros([m,m])
#PP_list=[]
#for i in range(m):
#    for j in range(m):
#        if rho[i,j] !=0:
#            PP_entry=P[i]*P[j]
#            PP_list.append(PP_entry)



#PP=np.zeros([m,m])
#for i in range(m):
#    for j in range(m):
#        PP[i,i]=PP_list[i]

##create g vector
#g_t=np.zeros([n+m,1])
#g_b=np.dot(PP,rho_vec)
#g=np.vstack((g_t,g_b))

###create solver to find steady states 
## fsolve complains about shape 
#def r(xyz):
#    x=xyz[:]
#    r=np.array(f+np.dot(A,x)+g)
#    return r
#x0=np.ones([8,1])

#xyz=fsolve(r,x0)



#implications of such a high tolerance?
#root = optimize.newton(r, 10*np.ones([8,1]),tol=100)

   



