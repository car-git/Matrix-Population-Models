# INSTALL PACKAGE "resample " IF NECESSARY:
#install.packages('resample')
library(resample)

#### PARAMETERS #####
nsim = 5000
c = 0
#alpha = 0.1
#Z = 1.6449
alpha = 0.05
Z=1.96
#alpha = 0.01
#Z = 2.5758
name_1 = 'LT3.csv' # <<<<< CHANGE OUTPUT FILE HERE 
#################




# FUNCTION FOR BOOTSTRAP
step1 = function(P,n,c,FF,S){
  muestreo = c()
  for (i in 1:S){
    muestreo = c(muestreo,rmultinom(1, n[i], P[i,]))
  }
  X = matrix(muestreo,nrow=S,byrow = T)
  Pnew = X/n
  U=t(Pnew[,-(S+1)])
  N=solve(diag(1,S)-U)
  lambda = Re(eigen(U+FF)$values[1])
  r = log(lambda)
  L = colSums(N)[1] # 1'Ne1
  R0 = colSums(FF%*%N)[1] # 1'FN,e1
  mu1 = colSums(FF%*%N%*%U%*%N)[1]/R0+c
  Tt = log(R0)/r
  x1=c(lambda,r,R0,mu1,Tt,L)
  return(x1)
}


file_readed = read.csv(name_1, header=TRUE) # <<<<< FUNCTION TO READ 

S = (ncol(file_readed)-2)/2
Etapas = colnames(file_readed)
Stages = Etapas[2:(2+S-1)]
A = data.matrix(file_readed)
time_units = A[,1]
Eggs = A[,(S+2)]
Nu = A[,2:(S+1)]
D = A[,(S+3):ncol(A)]

M = matrix(c(rep(0,length(time_units)*S)),ncol=S)

for (t in 1:(length(time_units)-1)){
  M[t+1,1] = Nu[t,1]-D[(t+1),1]-Nu[(t+1),1]
}

for (s in 2:S){
  for (t in 1:(length(time_units)-1)){
    M[t+1,s] = Nu[t,s]+M[(t+1),s-1]-D[(t+1),s]-Nu[(t+1),s]
  }
}

NN = colSums(Nu)
Dead = colSums(D)
Graduated = colSums(M)

G = Graduated/NN
R = Dead/NN
P = 1-G-R
U = diag(P)

for (i in 1:(S-1)){
  U[(i+1),i] = G[i]
}

fertility = sum(Eggs)/tail(NN,1)

FF = matrix(c(rep(0,length(U))),ncol=S)
FF[1,S] = fertility

# ************ STARTING CALCULATIONS ***********

###### EXPECTED VALUES
N=solve(diag(1,S)-U)
lambda = Re(eigen(U+FF)$values[1])
r = log(lambda)
L = colSums(N)[1] # 1'Ne1
R0 = colSums(FF%*%N)[1] # 1'FNe1
mu1 = colSums(FF%*%N%*%U%*%N)[1]/R0+c
Tt = log(R0)/r
Ex_vls = c(lambda,r,R0,mu1,Tt,L)

# Function to obtain confidence intervals
P1 = t(U)
P1 = cbind(P1,1-rowSums(P1))
muestra = c() 
for (i in 1:nsim){ ### Numbers of simulations
  muestra = c(muestra,step1(P1,NN,c,FF,S)) # Call the function
} 
Resultados = matrix(muestra,nrow=nsim,byrow = T) # Convert to a matrix
std_vls = sqrt(colVars(Resultados)) # Obtain standard deviation

# Life table parameters estimates

t = A[,1]
Tot = rowSums(Nu)[1] + rowSums(D)[1]
L = rowSums(A[,2:(S+1)])
L2 = L/Tot
M2 = Eggs/L

M2[is.nan(M2)] = 0

f = function(vx) t(exp(-vx*t)*L2)%*%M2-1
r_lt = uniroot(f, c(0,1))$root

la_lt = exp(r_lt)
Ro_lt = t(L2)%*%M2

mu_1_lt = t(t-1+c)%*%(L2*M2)/Ro_lt
T_lt = log(Ro_lt)/log(la_lt)
L_lt = sum(L2)

Lower = Ex_vls - Z*std_vls
Upper = Ex_vls + Z*std_vls

# Print the values
Table_5=matrix(c(NN,G,R,P,1/(1-P)),nrow=S)
rownames(Table_5) <- Stages
colnames(Table_5) <- c('Units', 'Gi', 'Ri', 'Pi', 'Ave')
as.table(Table_5)

Table_6 = matrix(c(la_lt,r_lt,Ro_lt,mu_1_lt,T_lt,L_lt,Ex_vls,std_vls),nrow=3,byrow=T) 
Table_6 = data.frame(Table_6)
colnames(Table_6) = c('λ','r','Ro','μ','T','L')
rownames(Table_6) = c('Life table','MPM','S.E')

Table_7 = matrix(c(Lower,Upper),nrow=2,byrow=T) 
Table_7 = data.frame(Table_7)
colnames(Table_7) = c('λ','r','Ro','μ','T','L')
rownames(Table_7) = c('Lower','Upper')

Table_5
Table_6
Table_7

