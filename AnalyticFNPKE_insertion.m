%Code that solves the Fractional Neutron Point Kinetics (FNPKE) for 
%step insertions.

%Authors: Cruz-López C.-A. (cacl.nucl@gmail.com)
%         Espinosa-Paredes G. (gepe@xanum.uam.mx)


format long g
%*******************************************************************
%************************** Input **********************************

%fractional order
alpha_f = 1;

%Vector with the standard decay lambda constants of the precursors
L=[0.0127 0.0317 0.115 0.311 1.4 3.87];

%Betas
Betas =[0.000285,0.0015975,0.00141,0.0030525,0.00096,0.000195];

%reactivity
rho = -0.0075

%Lambda_U
Lambda_U=0.0005;

% These lines do not need to be modified.
%Beta total, lambda^alpha, Lambda_U^alpha
L_f = L.^alpha_f;
LAM_f = Lambda_U^alpha_f;

%Initial conditions
n0 = 1;
C0 = (n0/LAM_f)*Betas./L_f

% The step variable is used for large times, with the purpose of
% avoiding possible numerical issues with the Mittag-Leffler. 
% For small values (of the order of 10 seconds), Target can be equal
% to the time step.

Target = 10
step = 10

Solution_n = [ ]
Solution_n=Insertion(Target,L_f,LAM_f,rho,Betas,step,alpha_f,C0);
filename = 'Output_insertion.xlsx';
xlswrite(filename,Solution_n)

function P1 = Poly_Coeff(L_f,LAM_f,rho,Betas)
C_P = [ ];
bet_tot = sum(Betas);
u = (rho-bet_tot)/LAM_f;
C_P(1:3)=[1 Suma(1,L_f)-u Suma(2,L_f)-u*Suma(1,L_f)-(1/LAM_f)*dot(L_f,Betas)];

for i=3:size(L_f,2)
    s1 = 0;
    for j=1:size(L_f,2)
        s1 = s1+L_f(j)*Betas(j)*Sumai(j,i-2,L_f);
    C_P(i+1)=Suma(i,L_f)-u*Suma(i-1,L_f)-(1/LAM_f)*s1;
    end
end
s2=0;
for k=1:size(L_f,2)
    s2 = s2+L_f(k)*Betas(k)*Sumai(k,size(L_f,2)-1,L_f);
end
C_P(size(L_f,2)+2)=-u*Suma(size(L_f,2),L_f)-(1/LAM_f)*s2;
P1 = C_P;
end

function P_d = Poly_Coeff_d(P)
exp = size(P,2)-1:-1:1;
P(1:size(P,2)-1);
P_d(1:size(P,2)-1)=P(1:size(P,2)-1).*exp;
end

function H1 = Poly_Coeff_H(L_f,LAM_f,rho,Betas,C0)
C_H = [ ];      
C_H(1)=dot(L_f,C0);
for j=2:size(L_f,2)
    s2=0;
    for i=1:size(L_f,2)
        s2 = s2+L_f(i)*C0(i)*Sumai(i,j-1,L_f);
    end
    C_H(j)=s2;
end
H1 = C_H;
end

function Q1 = Poly_Coeff_Q(L_f)
C_Q = [ ];
C_Q(1)=1;
for j=1:size(L_f,2)
    C_Q(j+1)=Suma(j,L_f);
end
Q1=C_Q;
end
    
%Sum with replacement
function S1i = Sumai(i,m,L)
COPY = L;
COPY(i)=[];
S1i=Suma(m,COPY);
end

%Sum with replacement
function S1 = Suma(m,L)
COMB = nchoosek(L,m);
Vec = prod(COMB,2);
S1 = sum(Vec);
end

function Sol_C = Solution_C(n0,P,P_d,Q,H,t,alpha_f,L_f,Betas,LAM_f,C0)
r_P = roots(P);
Solution_C = [ ];
for j=1:size(Betas,2)
    s1 = 0;
    for k=1:size(r_P,1)
        Num = n0*polyval(Q,r_P(k))+polyval(H,r_P(k));
        Quot = polyval(P_d,r_P(k))*(r_P(k)+L_f(j));
        diff = ml(r_P(k)*(t^(alpha_f)),alpha_f)-ml(-L_f(j)*(t^(alpha_f)),alpha_f);
        s1 = s1+Num*diff/Quot;
    end
    Solution_C(j)=Betas(j)*(s1/LAM_f)+C0(j)*ml(-L_f(j)*(t^(alpha_f)),alpha_f);

end
Sol_C = Solution_C;
end

function Sol = Solution(n0,P,P_d,Q,H,t,alpha_f)
r_P = roots(P);
dime = size(r_P,1);
s1 = 0;
for k=1:size(r_P,1)
    Num = n0*polyval(Q,r_P(k))+polyval(H,r_P(k));
    Quot = polyval(P_d,r_P(k));
    raiz = r_P(k);
    Mittag = ml(raiz*(t^alpha_f),alpha_f);
    
    s1 = s1+Num*Mittag/Quot;
end
Sol = s1;

end

function R =Insertion(Target,L_f,LAM_f,rho,Betas,h,alpha_f,C0)
Vec_sol = [ ]
bet_tot = sum(Betas);
t=h;
n0=1;
P = Poly_Coeff(L_f,LAM_f,rho,Betas);
Q = Poly_Coeff_Q(L_f);
P_d =  Poly_Coeff_d(P);

for m=1:Target/h
    solucion = 0;
    H = Poly_Coeff_H(L_f,LAM_f,rho,Betas,C0);
    n_f = Solution(n0,P,P_d,Q,H,t,alpha_f);
    C_f = Solution_C(n0,P,P_d,Q,H,t,alpha_f,L_f,Betas,LAM_f,C0);
    n0 = n_f;
    C0=C_f;
    m
    Vec_sol=[Vec_sol;t*m n_f C0(2)];
end
n0;
rho
R = Vec_sol
end
