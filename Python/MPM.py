 #!/usr/bin/env python3.6
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from tabulate import tabulate

# *********   Parameters   **********
fname = "LT3.csv"
nsim = 5000
c=0
# alpha = 0.1; Z=1.6449;
alpha = 0.05; Z=1.96;
# alpha=0.01; Z=2.5758
# ***********************************


def simul(U,v,u,I):
    N = np.linalg.inv(I-U)
    p1 = np.matmul(v,N)
    mu = np.matmul(p1,u).tolist()
    mu=mu[0][0]
    p1 = np.matmul(v,N)
    p2 = np.matmul(2*N-I,u)
    p3 = np.matmul(p1,p2)
    sigma = p3-mu*mu
    sigma = sigma[0][0]
    return mu,sigma

def estpars(U,F):
    S = np.shape(U)[0]
    I = np.eye(S)
    N = np.abs(np.linalg.inv(I-U))

    c=0
    u = np.array([np.ones(len(U))])
    e1 = np.array([np.zeros(len(U))])
    e1[0][0]=1
    e1=np.transpose(e1)
    FN = np.matmul(F,N)
    step1=np.matmul(u,FN)
    Roc=np.matmul(step1,e1)[0][0]


    UN = np.matmul(U,N);
    FNUN = np.matmul(FN,UN)
    step1 = np.matmul(u,FNUN)
    mu_1 = np.matmul(step1,e1)[0][0]/Roc + c

    L = np.sum(N[:,0])
    lam = np.max(np.real(np.linalg.eigvals(U+F)))
    r = np.log(lam)
    T = np.log(Roc)/np.log(lam)

    return lam,r,Roc,mu_1,T,L

def simula(U,F,n,nsim):

    experiments=[]

    Ut = np.transpose(U)
    for s in range(nsim):
        Uest = [] 
        row = 0 
        for i in Ut:
            p=list(i)
            p.append(1-sum(i))
            x = (np.random.multinomial(n[row], p, size=1)/n[row]).tolist()
            Uest.append(x[0][0:-1])
            row += 1
        Uest = np.transpose(Uest)
    
        try:
            [lam,r,Roc,mu_1,L,T] = estpars(Uest,F)
            experiments.append([lam,r,Roc,mu_1,L,T])
        except:
            pass
    
    # avsigma = np.mean(SIGMA)

    return experiments


def LifeTable(t,L,M,c,x_initial_guess):

    func = lambda x : np.sum(np.exp(-x*t)*L*M)-1
    r_LT = fsolve(func, x_initial_guess)[0]
    lam_LT = np.exp(r)
    Ro_LT = np.sum(L*M)
    mu_1_LT = np.sum((t-1+c)*L*M)/Ro_LT
    T_LT = np.log(Ro_LT)/np.log(lam_LT)
    L_LT = np.sum(L)
   

    return lam_LT,r_LT,Ro_LT,mu_1_LT,T_LT,L_LT



df = pd.read_csv(fname,header=0)
col_name=list(df.columns)
col_names = [x.replace('.','_') for x in col_name]

A = df.to_numpy()
S = int((A.shape[1]-2)/2)
Stages = col_names[1:S+1]
time_units = A[:,0]
Eggs=A[:,S+1]
Nu=A[:,1:S+1]
D=A[:,S+2::]
    
nrows= np.shape(A)[0]

M=np.zeros(np.shape(Nu))

for t in range(nrows-1):  
    M[t+1,0] = Nu[t,0] - D[t+1,0] - Nu[t+1,0]   

for s in range(1,S):
    for t in range(nrows-1):     
        M[t+1,s] = Nu[t,s] + M[t+1,s-1] - D[t+1,s] - Nu[t+1,s]

    # wait = input("Press Enter to continue.")
NN = Nu.sum(axis=0)

Dead = D.sum(axis=0)
Graduated = M.sum(axis=0)

R = Dead/NN
G = Graduated/NN
P = 1-R-G
U=np.diag(P)
for i in range(S-1):
    U[i+1,i] = G[i]

fertility = np.sum(A[:,S+1])/NN[-1]
F=np.zeros(np.shape(U))
F[0,-1]=fertility
[lam,r,Roc,mu_1,T,L] = estpars(U,F)


expe = simula(U,F,NN,nsim)

med=np.mean(expe,axis=0)
SE=np.std(expe,axis=0)


# Life table calculations

N = np.sum(Nu[0])+np.sum(D[0])

L_LT=np.sum(Nu,axis=1)
M_LT = Eggs
pos=np.where(L_LT>0)

L_LT = L_LT[pos]
M_LT = M_LT[pos]
M_LT = M_LT/L_LT
L_LT = L_LT/N
t_LT = time_units[pos]
ini_guess = r
res_LT = LifeTable(t_LT,L_LT,M_LT,c,ini_guess)


# Preparing tables to print
I = np.eye(S)
N = np.abs(np.linalg.inv(I-U))

cols5 = ['Units', 'Gi', 'Ri', 'Pi', 'Ave']
Table_5=[]
for i in range(S):
    Table_5.append([Stages[i],NN[i],np.round(G[i],4),np.round(P[i],4),np.round(1/(1-P[i]),4)])

cols6 = ['Method', 'lambda', 'r', 'Ro', 'mu_1','T','L']
Table_6=[]
a=["Life table", np.round(res_LT[0],4),np.round(res_LT[1],4),np.round(res_LT[2],4),np.round(res_LT[3],4),np.round(res_LT[4],4),np.round(res_LT[5],4)]
Table_6.append(a)
b=["MPM",np.round(lam,4),np.round(r,4),np.round(Roc,4),np.round(mu_1,4),np.round(T,4),np.round(L,4)]

Table_6.append(b)

cols7 = ['', 'lambda', 'r', 'Ro', 'mu_1','T','L']
Table_7=[]
low = np.round(med-1.96*SE,4)
upp = np.round(med+1.96*SE,4)

a = ["Lower",np.round(low[0],4),np.round(low[1],4),np.round(low[2],4),np.round(low[3],4),np.round(low[4],4),np.round(low[5],4)]

Table_7.append(a)
b = ["Upper",np.round(upp[0],4),np.round(upp[1],4),np.round(upp[2],4),np.round(upp[3],4),np.round(upp[4],4),np.round(upp[5],4)]

Table_7.append(b)

print("  ")
print("  ")
print("  ")
print("                         RESULTS")
print("  ")
print("  ")
print('Matrix U : ')
print("  ")
print(np.round(U,4))
print("  ")
print("  ")
print('Matrix F : ')
print("  ")
print(np.round(F,4))
print("  ")
print("  ")
print('Matrix N : ')
print("  ")
print(np.round(N,4))
print("  ")
print("  ")
print("Number of simulations for boostrap: ")
print("  ")
print(nsim)
print("  ")
print("  ")
print("Table_5: ")
print("  ")
print(tabulate(Table_5, headers=cols5))
print("  ")
print("  ")
print("Table_6: ")
print("  ")
print(tabulate(Table_6, headers=cols6))
print("  ")
print("  ")
print("Table_7: ")
print("  ")
print(tabulate(Table_7, headers=cols7))


df = pd.DataFrame (Table_5, columns=cols5)
df.to_csv("Table_5.csv",index=False)

df = pd.DataFrame (Table_6, columns=cols6)
df.to_csv("Table_6.csv",index=False)

df = pd.DataFrame (Table_7, columns=cols7)
df.to_csv("Table_7.csv",index=False)




