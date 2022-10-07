#!/usr/bin/env wolframscript
(* ::Package:: *)

fname="LT3.csv";
nsim=5000; (* 5,000 simulations is usually enough *)
c=0;   (* c is usually set to zero *)

(* Select the alpha value for CI (default is 0.05) *)

(*alpha = 0.1; Z=1.6449;*)
alpha = 0.05; Z=1.96;
(*alpha=0.01; Z=2.5758*)

        (* Do not make changes below this line *)

SetDirectory@NotebookDirectory[];
data=Import[fname,"Data","HeaderLines"->1];
data=Import[fname,"Data"];
S=((Dimensions[data[[1]]]-2)/2)[[1]];
St=data[[1]];
Stages=Take[St,{2,S+1}];
nrows=Dimensions[data][[1]];
A=Take[data,{2,nrows}];
nrows = nrows-1;
TimeUnits=A[[All,1]];
Eggs=A[[All,S+2]];
Nu=A[[All,2;;S+1]];
De=A[[All,S+3;;]];
M=ConstantArray[0,{nrows,S}];

For[t=1,t<nrows,t++,
M[[t+1,1]]=Nu[[t,1]]-De[[t+1,1]]-Nu[[t+1,1]]
]

For[s=2,s<=S,s++,
For[t=1,t<nrows,t++,
M[[t+1,s]]=Nu[[t,s]]+M[[t+1,s-1]]-De[[t+1,s]]-Nu[[t+1,s]]
]]

Nt=Total[Nu];
Dead=Total[De];
Graduated=Total[M];
G=Graduated/Nt ;
R=Dead/Nt ;
P=1-G-R;
U=DiagonalMatrix[P];
For[s=1,s<S,s++,
U[[s+1,s]]=G[[s]]
]
fertility = Total[Eggs]/Nt[[S]];
F=ConstantArray[0,{S,S}];
F[[1,S]]=fertility;
Id=IdentityMatrix[S];
NN=Inverse[Id-U];
ET = Diagonal[NN];
e1=ConstantArray[0,{S,1}];
 e1[[1,1]]=1;
u=ConstantArray[1,{1,S}];
Ro=u.F.NN.e1;
Ro = Ro[[1]][[1]];
mu1=u.F.NN.U.NN.e1/Ro+c;
mu1 = mu1[[1]][[1]];
L = Total[NN[[All,1]]];
lam=Eigenvalues[U+F];
lam=Select[lam,Im[#]==0&];
\[Lambda]=Max[lam]//N;
r=Log[\[Lambda]];
T=Log[Ro]/Log[\[Lambda]];
f=1-Total[Transpose[U],{2}];
P=MapThread[Append,{Transpose[U],f}];

res=ConstantArray[0,{nsim,6}];
For[i=1,i<nsim+1,i++,
{Pest=ConstantArray[0,{S,S}],
For[j=1,j<S+1,j++,
{Y = RandomVariate[MultinomialDistribution[Nt[[j]],P[[j]]],1],
Pest[[j]]=Y[[1]]/Nt[[j]]}
],
Ut=Transpose[Pest],
Ut=Ut[[1;;S,All]],
NNt=Inverse[Id-Ut];
ETt = Diagonal[NNt];
Rot=u.F.NNt.e1;
Rot= Rot[[1]][[1]]//N;
mu1t=u.F.NNt.Ut.NNt.e1/Rot+c;
mu1t = mu1t[[1]][[1]];
Lt = Total[NNt[[All,1]]]//N;
lamt=Eigenvalues[Ut+F];
lamt=Select[lamt,Im[#]==0&];
\[Lambda]t=Max[lamt]//N;
rt=Log[\[Lambda]t];
Tt=Log[Rot]/Log[\[Lambda]t];
v={\[Lambda]t,rt,Rot,mu1t,Tt,Lt};
res[[i]]=N[v];
}
]
me=Mean[res];
sd=StandardDeviation[res];

Ntt=Total[Nu[[1,All]]]+Total[De[[1,All]]];
L=Total[Nu,{2}];
pos=Position[L,0];
If[Length[pos]>0,pos=Min[pos],pos=0];
If[pos>0,L=Take[L,{1,pos-1}]];
If[pos>0,Eggs=Take[Eggs,{1,pos-1}]];
If[pos>0,TimeUnits=Take[TimeUnits,{1,pos-1}]];
m=Eggs/L;
L2=L/Ntt;
sol=FindRoot[Total[Exp[-r2*TimeUnits]*L2*m]-1,{r2,r}];
r2=r2/.sol;
\[Lambda]2 = Exp[r2];
Ro2=L2.m //N;
\[Mu]2 = Total[(TimeUnits-1+c)*L2*m]/Ro2;
T2=Log[Ro2]/Log[\[Lambda]2];
LL2=Total[L2]//N;
low=me-Z*sd;
upp=me+Z*sd;

Matrix U :
MatrixForm[U]//N

Matrix F:
MatrixForm[F]//N

Matrix N:
MatrixForm[NN]//N





Table 5 :
G=N[G];
R=N[R];P=1-G-R;T=1/(1-P);
Print["-----------------------------------------------------------------------"];
Print["Stage Units           Gi             Ri             Pi            Ti"];
Print["-----------------------------------------------------------------------"];
For[i=1,i<S+1,i++,
{Print[i,"       ",N[Nt[[i]],4],"       ",N[G[[i]],4],"       ",N[R[[i]],4],"       ",N[P[[i]],4],"       ",N[T[[i]],4]],
}
]
Print["-----------------------------------------------------------------------"];





Table 6 :
Print["--------------------------------------------------------------------------------"];
Print["Method           \[Lambda]             r         Ro          \[Mu]1          T          L"];
Print["--------------------------------------------------------------------------------"];
Print["Life Table    ",N[\[Lambda]2],"     ",r2,"     ",Ro2,"     ",\[Mu]2,"     ",T2,"     ",LL2 ];
Print["MPM           ",\[Lambda]t,"     ",rt,"     ",Rot,"     ",mu1t,"     ",Tt,"     ",Lt ];
Print["--------------------------------------------------------------------------------"];



Table 7 :
Print["--------------------------------------------------------------------------------"];
Print["              \[Lambda]             r         Ro          \[Mu]1          T          L"];
Print["--------------------------------------------------------------------------------"];
Print["Lower      ",low[[1]],"     ",low[[2]],"     ",low[[3]],"     ",low[[4]],"     ",low[[5]],"     ",low[[6]] ];
Print["Upper      ",upp[[1]],"     ",upp[[2]],"     ",upp[[3]],"     ",upp[[4]],"     ",upp[[5]],"     ",upp[[6]] ];
Print["--------------------------------------------------------------------------------"];



