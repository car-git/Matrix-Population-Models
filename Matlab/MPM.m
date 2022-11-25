function MPM(datafile)

% ﻿*********   Parameters   **********
nsim = 5000;
c=0;
%alpha = 0.1; Z=1.6449;
alpha = 0.05; Z=1.96;
%alpha=0.01; Z=2.5758
% *****************************


T = readtable(datafile);
S = (size(T,2)-2)/2;
Stages = T.Properties.VariableNames(2:2+S-1)';
A=table2array(T);

time_units = A(:,1);
Eggs=A(:,S+2);
Nu=A(:,2:S+1);
D=A(:,S+3:end);

% First stage movements
for t = 1 : size(time_units,1)-1
    M(t+1,1)=Nu(t,1) - D(t+1,1) - Nu(t+1,1);
end


% Rest of stage movements
for s=2:S
    for t = 1 : size(time_units,1)-1
        M(t+1,s)=Nu(t,s) + M(t+1,s-1) - D(t+1,s) - Nu(t+1,s);
    end
end

% Now begins analysis per stage

NN = sum(Nu);
Dead = sum(D);
Graduated = sum(M);

G=Graduated./NN;
R = Dead./NN;
P = 1-G-R;

U=diag(P);

for i=1:S-1
    U(i+1,i)=G(i);
end
clc
fprintf('                   Transition matrix U:  \n') 
U

Table_5=table();
Table_5.Units=NN';
Table_5.Gi = round(G,4)';
Table_5.Ri = round(R,4)';
Table_5.Pi = round(P,4)';
Table_5.Ave = round(1./(1-P),4)';
Table_5.Properties.RowNames=Stages;


% Calculation of LHT (MPM)

%fertility = sum(Eggs)/sum(Nu(:,end));
fertility = sum(Eggs)/NN(end)


fprintf('                   Fertility matrix F:  \n') 
F=zeros(size(U));
F(1,end)=fertility
I=eye(S);
fprintf('                   Fundamental matrix N:  \n')

N=inv(I-U)
ET=diag(N);
%fprintf('                  Average residence times (diagonal of fundamental matrix):  \n')



e1=zeros(S,1); e1(1,1)=1;
u=ones(S,1);
Roc=u'*F*N*e1;
muc=(1/Roc)*u'*F*N*U*N*e1+c;
Lc=N(:,1)';
lambdac=max(eig(U+F));
T=log(Roc)/log(lambdac);


Ro=sum(Eggs)/Nu(1,1);

y=[lambdac log(lambdac) Roc muc T  sum(Lc)];

sd = CIn(U,F,NN,nsim);

out=[y;sd] ;
fprintf('Number of simulations for boostrap:  \n')
nsim

Lo=out(1,:)-Z*out(2,:);
Up=out(1,:)+Z*out(2,:);
ci=[Lo;Up];


% Life table parameter estimates

Tot=sum(Nu(1,:))+sum(D(1,:));

t=A(:,1);
L=sum(A(:,2:S+1),2);
L2=L/Tot;
M2=Eggs./L;

M2(isnan(M2))=0;

save data t L2 M2


%Define function handle
f = @(x) lambd(x,t,L2,M2);
% initial value is lambdac
r_lt = fzero(f,log(lambdac))


la_lt = exp(r_lt);
Ro_lt = L2'*M2;

mu_1_lt = (t-1+c)'*(L2.*M2)/Ro_lt;
T_lt=log(Ro_lt)/log(la_lt);
L_lt = sum(L2);


Table_6=table();
Table_6.lambda=round([la_lt lambdac sd(1)]',3);
Table_6.Properties.VariableNames(1) = "λ";
Table_6.r=round([r_lt log(lambdac) sd(2)]',3);
Table_6.Ro=round([Ro_lt Roc sd(3)]',3);
Table_6.muc=round([mu_1_lt muc sd(4)]',3);
Table_6.Properties.VariableNames(4) = "μ";
Table_6.T=round([T_lt T sd(5)]',3);
Table_6.L=round([L_lt sum(Lc) sd(6)]',3);
Table_6.Properties.RowNames={'Life Table','MPM','S.E.'};



Table_7=table();
Table_7.lambda=round(ci(:,1),3);
Table_7.Properties.VariableNames(1) = "λ";
Table_7.r=round(ci(:,2),3);
Table_7.Ro=round(ci(:,3),3);
Table_7.muc=round(ci(:,4),3);
Table_7.Properties.VariableNames(4) = "μ";
Table_7.T=round(ci(:,5),3);
Table_7.L=round(ci(:,6),3);
Table_7.Properties.RowNames={'Lower','Upper'};


Table_5
Table_6
Table_7
save('results','U','F','N','Table_5','Table_6','Table_7')



function out = lambd(x,t,L,M)
    out = sum(exp(-x*t).*L.*M)-1;


function Res = CIn(U,F,n,nsim)

S = size(U,1);
u=ones(S,1);
f = u-sum(U',2);
P = [U' f];

for i=1: nsim
    Pest=zeros(S+1);
    Y = zeros(S+1);

    Y=mnrnd(n',P);
    for j=1:S
        Pest(j,:)=Y(j,:)/n(j);
    end
    Ui=Pest(1:S,1:S);
    Ui=Ui';
    c=0; % (Using c=0)
    e1=zeros(S,1); e1(1,1)=1;

    I=eye(S);
    N = inv(I-Ui);
    R=u'*F*N*e1;
    mu=(1/R)*u'*F*N*Ui*N*e1+c;
    L=sum(N(:,1));
    lambda=max(eig(Ui+F));
    T=log(R)/log(lambda);
    r=log(lambda);

    W(i,:) = [lambda r R mu T L];

end
Res = std(W);




