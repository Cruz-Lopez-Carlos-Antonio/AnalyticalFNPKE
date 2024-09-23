# Analytical Solution of the Fractional Neutron Point Kinetic Equations. 
The present repository contains the MATLAB codes developed to solve the Fractional Neutron Point Kinetic Equations (FNPKE). These codes were described in the article *A New Solution of the Fractional Neutron Point Kinetics Equations using Symmetry and the Heaviside’s expansion formula*, which was recently submitted to the **Progress in Nuclear Energy** journal. 

The programs are licensed under a Creative Commons Attribution 4.0 International License: http://creativecommons.org/licenses/by/4.0/

Authors: Carlos-Antonio Cruz-López (cacl.nucl@gmail.com), Gilberto Espinosa-Paredes (gepe@xanum.uam.mx)

Mathematical and algorithmical generalities of the codes are described in the following lines with the purpose to provide some insight of the developed work. Nevertheless, a more detailed discussion is provided in the submitted article.
## Software specifications and requirements. 
The AnalyticFNPKE codes were written in the MATLAB programming language R2021a, but older versions can be used, from the R2012b one onward. The reported examples and results were obtained in a 3.8 GHz desktop computer, with 32 Gb of RAM and under a Windows 11 environment. The developed codes require the following script, which was written by Roberto Garrappa (2015, 2024):
- [x] ml.m

This last file can be freely download from the mathworks site in the following [link](https://www.mathworks.com/matlabcentral/fileexchange/48154-the-mittag-leffler-function).

## Financial Support.
The authors appreciate the financial support received from the Consejo Nacional de Humanidades, Ciencia y Tecnología (CONAHCYT), under the program *Estancias Posdoctorales por México, 2022*, with the project entitled: *Desarrollo de modelos fenomenológicos energéticos de orden fraccional, para la optimización y simulación en reactores nucleares de potencia*, by which the present development was possible.

## Index of the Repository
1. [Mathematical description of the problem.](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticalFNPKE#1-mathematical-description-of-the-problem)
1. [Laplace transform of the system.](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticalFNPKE#2-laplace-transform-of-the-system)
1. [Analytical solution and Polynomials.](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticalFNPKE#3-analytical-solution-and-polynomials)
1. [Equivalence between the integer's solution and the fractional one.](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticalFNPKE/blob/main/README.md#4equivalence-between-the-integers-solution-and-the-fractional-one)
   - [4.1 Equivalence between the sums.](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticalFNPKE/blob/main/README.md#41-equivalence-between-the-sums)
   - [ 4.2 Equivalence between the Polynomials.](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticalFNPKE/blob/main/README.md#42-equivalence-between-the-polynomials)
   - [ 4.3 Equivalence between the Mittag-Leffler function and the exponential one.](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticalFNPKE#43-equivalence-between-the-mittag-leffler-function-and-the-exponential-one)
3. [Algorithmical Implementation.](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticalFNPKE#5-algorithmical-implementation)
   - [5.1 Requeriments.](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticalFNPKE#51-requirements)
   - [5.2 From Python to MATLAB.](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticalFNPKE#52-from-python-to-matlab)
   - [5.3 Sums.](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticalFNPKE#53-sums)
   - [5.4 Shifted sums.](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticalFNPKE#54-shifted-sums)
   - [5.5 Polynomials.](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticalFNPKE#55-polynomials)
1. [AnalyticFNPKE_Insertion.m](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticalFNPKE#6-analyticfnpke_insertionm)
   - [6.1 Verification of the AnalyticFNPKE_Insertion.m](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticalFNPKE#61-verification-of-the-analyticfnpke-insertion.m)
   - [6.2 Example of a calculation with the AnalyticFNPKE_Insertion.m](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticalFNPKE#62-calculations-for-the-fractional-order-alpha--09)
1. [AnalyticFNPKE_Ramp.m](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticalFNPKE#7-analyticfnpke_rampm)
   

## 1. Mathematical description of the problem
Following the ideas of Ray and Patra (2014), a compartmental version of the The Fractional Neutron Point Kinetic Equations (FNPKE), with $K$ groups of precursors of delayed neutrons, can be written as follows:
$$D_C^\alpha n\left(t\right)=\frac{\rho\left(t\right)-\beta}{\Lambda^\alpha}n\left(t\right)+\sum_{i=1}^{K}{\lambda_i^\alpha C_i(t)}, \tag{1}$$

$$D_C^\alpha C\left(t\right)=\frac{\beta_i}{\Lambda^\alpha}n\left(t\right)-\lambda_i^\alpha C_i\left(t\right),\ \ \ \ \ 1\le\ i\le\ K, \tag{2}$$

where $n(t)$ denotes the neutron density and $\rho(t)$  the reactivity. $C(t),\lambda_k,\beta_k$ are the concentration, the decay constant and the fraction of the $k$-group of precursors of the delayed neutrons, respectively. $\Lambda$ represents the prompt neutron generation time,  $\beta$ is given by:

$$\sum_{k=1}^{K}\beta_k=\beta, \tag{3}$$

and the operator $D_C^\alpha$ denotes the fractional derivative of Caputo, defined as:
$$D_C^\alpha f\left(t\right)=\frac{1}{\Gamma(n-\alpha)}\int_{0}^{t}{f^{\left(n\right)}\left(\tau\right)\left(t-\tau\right)^{n-\alpha-1}d\tau}, \tag{4}$$

where $\alpha$ is the fractional order, $\Gamma(\cdot)$ is the Gamma function, and $n \in \mathbb{N}$ fulfills the following inequality:

$$n-1 \leq \alpha \< n. \tag{5}$$

This fractional version is different from the one proposed by Espinosa-Paredes et al. (2011), because it was developed in terms of a fractional mass balance approach, instead of a transport theory's one. The physical implications are discussed in detail in the submitted paper, but it can be understood, essentially, as a fractional compartmental model similar to the one developed in the Pharmacokinetic field (Dokoumetzidis, 2010).
Nahla developed an analytic solution of a system related to Eq. (1) and Eq.(2) (2017) using the Laplace transform and a matrix approach. We proposed a different approach using algebraic theory of equations and a procedure that uses symmetry. 

## 2. Laplace transform of the system.
The Laplace transform of the Caputo's derivative fulfills the following relationship (Ishteva, 2005, p. 21):

$$\mathcal{L} \\{D_C^\alpha (t),s \\}=s^\alpha F\left(s\right)-\sum_{k=0}^{n-1}{s^{\alpha-k-1}f^{\left(k\right)}(0)}, \tag{6}$$

where $F(s)=\mathcal{L} \\{f(t),s \\}$ denotes the standard Laplace transform of $f(t)$ and $f^{k}(0)$ are the initial conditions of the function. Applying the last relationship on both sides of Eq. (1) and Eq. (2), the following equations are obtained:

$$s^\alpha\widetilde{n}\left(s\right)-s^{\alpha-1}n\left(0\right)=\frac{\rho-\beta}{\Lambda^\alpha}\widetilde{n}\left(t\right)+\sum_{i=1}^{K}{\lambda_i^\alpha{\widetilde{C}}_i}\left(s\right), \tag{7}$$

$$s^\alpha{\widetilde{C}}_i\left(s\right)-s^{\alpha-1}C_i\left(0\right)=\frac{\beta_i}{\Lambda^\alpha}\widetilde{n}\left(s\right)-\lambda_i^\alpha{\widetilde{C}}_i(s). \tag{8}$$

Using a new procedure based on the symmetry between this fractional system and the integer one, it is possible to rewrite the last systems of equation as follows:

$$\widetilde{n}\left(s\right)=\sum_{j=1}^{K+1}{A_j\frac{s^{\alpha-1}}{s^\alpha-p_{j,f}}}, \tag{9}$$

$${\widetilde{C}}_{i}(s) =\frac{\beta_i}{\Lambda^\alpha} \sum _{j=1}^{K+1}A_j \frac{s^{\alpha-1}}{s^\alpha-p _{j,f}} \frac{1}{s^\alpha+\lambda_i^\alpha} +C_i\left(0\right)\frac{s^{\alpha-1}}{s^\alpha+\lambda_i^\alpha}, \tag{10}$$

where the coefficients $A_j$ are constants defined as:
$$A_j=\frac{n\left(0\right)Q_f\left(p_{j,f}\right)+H_f(p_{j,f})}{P_f^\prime(p_{j,f})}, \tag{11}$$

being $Q_f(s), H_f(s)$ and $P'_ f(s)=dP_f(s)/ds$ polynomials that are evaluated at the set of real numbers $p_{j,f}$, which in turn are the roots of the $P_f(s)$ polynomial. **A more detailed explanation about the deduction of the last expressions and the way in which the quotient of polynomials was obtained, is provided in the submitted paper.** Only the main mathematical aspects are discussed in the present repository, with the purpose to develop a computational algorithm.

## 3. Analytical solution and Polynomials.

Using Eq. (9) and Eq. (10), it is possible to find the solutions of $n(t)$ and $C_i(0)$ as follows:

$$n\left(t\right)=\sum_{j=1}^{K+1}A_jE_{\alpha,1}(p_{j,f}t^\alpha), \tag{12}$$

$$C_i\left(t\right)=\frac{\beta_i}{\Lambda^\alpha}\sum_{j=1}^{K+1}A_j\frac{E_{\alpha,1}\left(p_{j,f}t^\alpha\right)-E_{\alpha,1}(-\lambda_i^\alpha t^\alpha)}{p_{j,f}+\lambda_i^\alpha}+C_i\left(0\right)E_{\alpha,1}(-\lambda_i^\alpha t^\alpha), \tag{13}$$

where $E_{\alpha,1}(z)$ is the Mittag-Leffler function, which can be defined as (Gorenflo et al., 2020, p. 64):

$$E_{\alpha,\beta}\left(z\right)=\sum_{k=0}^{\infty}\frac{z^k}{\Gamma(k\alpha+\beta)},\ \tag{14}$$ 

with $\mathfrak{R}\left(\alpha\right)>0,\ \beta\in\mathbb{C}$. On the other hand, explicit expressions of the Polynomials related to coefficients $A_j$ are given by:
$$P_f\left(s\right)=s^{K+1}+\left(S_{1,K}-u\right)s^K+\left(S_{2,K}-uS_{1,K}-\frac{1}{\Lambda^\alpha}\sum_{i=1}^{K}{\lambda_i^\alpha\beta_i}\right)s^{K-1}$$

$$+\sum_{i=3}^{K}{\left(S_{i,K}-uS_{i-1,K}-\frac{1}{\Lambda^\alpha}\sum_{j=1}^{K}{\lambda_j^\alpha\beta_jS_{i-2,K-1}^j}\right)s^{K+1-i}-uS_{K,K}}$$

$$-\frac{1}{\Lambda^\alpha}\sum_{i=1}^{K}{\lambda_i^\alpha\beta_iS_{K-1,K-1}^i}; \tag{15}$$

$$H_f(s)=\sum_{i=1}^{K}\lambda_i^\alpha C_i(0)s^{K-1}+\sum_{j=2}^{K} \sum_{i=1}^{K}{\lambda_i^\alpha C_i(0) S_{j-1,K-1}^i}s^{K-j}; \tag{16}$$ 

and:
$$Q_f\left(s\right)=s^K+\sum_{j=1}^{K}{S_{j,K}s^{K-j}}, \tag{17}$$

where the sums $S_{m,n}$ and $S_{m,n}^i$ are defined as:

$$ S_{m,n}=\sum_{k_1=1}^{n-m+1}\sum_{k_2=k_1+1}^{n-m+2} \cdots \sum_{k_m=k_{m-1}+1}^{n} {\lambda_{k_1}^{\alpha} \lambda_{k_2}^{\alpha} \cdots\lambda_{k_m}^{\alpha}}, \tag{18}$$

$$S_{m,n}^i=\sum_{k_1=1,\ k_1\neq i}^{n-m+1}{\ \sum_{k_2=k_1+1,k_2\neq i}^{n-m+2}\cdots}\sum_{k_m=k_{m-1}+1,\ k_m\neq i}^{n}{\lambda_{k_1}^{\alpha} \lambda_{k_2}^{\alpha} \cdots\lambda_{k_m}^{\alpha}}, \tag{19}$$

and:
$$u=\frac{\rho-\beta}{\Lambda^\alpha} \tag{20}$$


## 4.Equivalence between the integer's solution and the fractional one. 

It is possible to show that the developed solution can be reduced to the expressions related to the integer case. These last were developed in the work *A New Simplified Analytical Solution to Solve the Neutron Point Kinetics Equations Using the Laplace Transform Method*, which was publised in the journal *Computer Physics Communications*. We made a particular repository for the integer case that can be consulted in the following [link](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE/tree/main).

### 4.1 Equivalence between the sums.
It is possible to understand the sums that were defined in Eq. (18) and Eq. (19) as functions over a set. In fact, it is necessary to define:


$$\mathcal{S}_m:\Omega\rightarrow\mathbb{R}, m\in \mathbb{N}, m\leq n,$$

where $\Omega$ is a set of real numbers given by $\Omega = \set{a_1,a_2,\cdots,a_n}$, and:

$$S_m(\Omega)=S_m(\set{a_1,a_2,\ldots,a_n})=\sum_{k_1=1}^{n-m+1} \sum_{k_2=k_1+1}^{n-m+2}\cdots \sum_{k_m=k_{m-1}+1}^{n}{a_{k_1}a_{k_2}\cdots a_{k_m}}. \tag{21}$$

Therefore, the sums can be rewritten as follows:
 
$$S_m(\set{\lambda_1^\alpha,\lambda_2^\alpha,\ldots,\lambda_n^\alpha})=S_m(\set{\lambda_i^\alpha\|1\leq i \leq\ n})\ = S_{m,n}. \tag{22}$$

Even more, using this function defined over sets, it is possible to compute the sum given in Eq. (19), that has a restriction over the indexes, as follows:
$$S_{m,n}^i=S_m(\Omega^i)=S_m (\set{\lambda_1^\alpha,\lambda_2^\alpha,\ldots,\lambda_{i-1}^\alpha,\lambda_{i+1}^\alpha,\ldots,\lambda_n^\alpha})$$
$$=S_m(\set{\lambda_k^\alpha |1\leq k\leq n,\ k\neq i}). \tag{23}$$

Using this notation, it follows that:

$$\lim_{\alpha \rightarrow 1} \underbrace{S_m(\Omega)}_ {\mathrm{Fractional\ case}}= S_m(\lim_{\alpha \rightarrow 1}{\Omega})=S_m(\lim_{\alpha \rightarrow 1} \set{\lambda_1^\alpha,\lambda_2^\alpha,\cdots,\lambda_n^\alpha})$$

$$= S_m(\set{\lambda_1,\lambda_2,\cdots,\lambda_n}) =\underbrace{S_m(\set{\lambda_1,\lambda_2,\ldots,\lambda_n}}_{\mathrm{Integer\ case}}. \tag{24}$$

This last expression is identical to one that is given for the integer case, as it can be corroborated in the following [link](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE#4-simplification-of-the-polynomials). A similar proof can be done for the case of $S_{m,n}^i$. 
### 4.2 Equivalence between the Polynomials. 

Once the equivalcente between the sums has been proved, it is possible to show the equivalcence of the Polynomials. This last can be done in a straighforward way observing that the polynomials can be written in terms of two variables as follows:

$$\lim_{\alpha\rightarrow1}{P_f(s,\alpha)}=\lim_{\alpha\rightarrow1}({s^{K+1}+\left(\fbox{$S_{1,K}$}-\frac{\rho-\beta}{\Lambda^{\fbox{$\alpha$}}}\right)s^K}$$

$$+ \left(\fbox{$S_{2,K}$}-\frac{\rho-\beta}{\Lambda^{\fbox{$\alpha$}}}\fbox{$S_{1,K}$}-\frac{1}{\Lambda^{\fbox{$\alpha$}}}\sum_{i=1}^{K}{\lambda_i^{\fbox{$\alpha$}}\beta_i}\right)s^{K-1}$$

$$+\sum_{i=3}^{K}{\left(\fbox{$S_{i,K}$}-\frac{\rho-\beta}{\Lambda^{\fbox{$\alpha$}}}\fbox{$S_{i-1,K}$}-\frac{1}{\Lambda^{\fbox{$\alpha$}}}\sum_{j=1}^{K}{\lambda_j^{\fbox{$\alpha$}}\beta_j\fbox{$S_{i-2,K-1}^j$}}\right)s^{K+1-i}-\frac{\rho-\beta}{\Lambda^{\fbox{$\alpha$}}}\fbox{$S_{K,K}$}}$$

$$-\frac{\rho-\beta}{\Lambda^{\fbox{$\alpha$}}}\fbox{$S_{K,K}$}-\frac{1}{\Lambda^{\fbox{$\alpha$}}}\sum_{i=1}^{K}{\lambda_i^{\fbox{$\alpha$}}\beta_i\fbox{$S_{K-1,K-1}^i$}}. \tag{25}$$

All the parts that depend on the fractional order $\alpha$ have been enclosed in a box, with the purpose to show the explicit dependence of that parameter. From Eq. (25) and Eq. (25) |it follows that:

$$\lim_{\alpha \rightarrow 1} \underbrace{P_f (s,\alpha)} _ {\mathrm{Fractional\ case}}=\underbrace{P(s)}_{\mathrm{Integer\ case}}, \tag{26}$$

where the form of the polynomial $P(s)$ can be consulted in the following [link](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE#4-simplification-of-the-polynomials). Similar proofs can be carried out for the polynomials $H(s)$ and $Q(s)$. 
### 4.3 Equivalence between the Mittag-Leffler function and the exponential one. 

From the definition given in Eq. (14), it follows that:

$$\lim_{\alpha \rightarrow 1} E_{\alpha,1}(z)=\sum_{k=0}^{\infty} \frac{z^k}{\Gamma(\alpha k+1)}=\sum_{k=0}^{\infty} \lim_{\alpha \rightarrow 1} \frac{z^k}{\Gamma(\alpha k+1)} \tag{27}$$

$$=\sum_{k=0}^{\infty} \frac{z^k}{\Gamma(k+1)}=\sum_{k=0}^{\infty}\frac{z^k}{k!}=\exp(z)$$

where the limit has been interchanged with the series due to the absolute convergence of the Mittag-Leffler function. Combining the results given in **Section 4.1**, **Section 4.2** and in the present section, it follows that the analytic solution of the FNPKE is obtained from the analytical solution of the fractional version. This result is relevant in terms of the algorithmical implementation, because allows extending all the codes that were developed for the NPKE (which can be consulted [here](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE#4-simplification-of-the-polynomials).

## 5. Algorithmical implementation. 
### 5.1 Requirements.

The present codes require, for their execution, of the **ml.m** code developed by Garrappa (2015, 2024), which can be freely download in the following [link](https://www.mathworks.com/matlabcentral/fileexchange/48154-the-mittag-leffler-function). It is worth mentioning that the use of such code is subject to a licence that must to be consulted. 
In order to use such code it is necessary to include it in the folder where the AnalyticFNPKE codes are saved. In the following image such procedure is showed:

![image](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticalFNPKE/assets/139827225/c370713f-75a2-4a4a-9bba-e0dfb8885cfc)
### 5.2 From Python to MATLAB.

Due to the symmetry and equivalence between the fractional and the integer case (when $\alpha =1$), in principle the codes that were developed for the NPKE in the Python language can be used, with a slightly modifications, to solve the FNPKE. Unfortunately there is an important obstacle that prevents that: the computation of the Mittag-Leffler function. 
This function requires for special algorithms to be accurately computed, because its standard definition given in Eq. (14) in terms of power series has a very slow convergence, requiring several terms to have an adequate result.
The main problem is related to the following term:
$$\Gamma(k\alpha+\beta) \tag{28}$$

which exhibits numerical issues for greater values of $k$, as Ortigueira et al. (2019) has pointed out. Therefore, it is necessary to use a more advanced method to compute such function, as the one provided by Roberto Garrappa (2015, 2024). Nevertheless, this author implements its algorithm in the MATLAB programming language, being a non-trivial task to migrate it to Python. Therefore it is necessary to change of programming language in order to include it.

### 5.3 Sums.
The first step of the algorithmical implementation oconsists of computing the sum given in the Eq. (18) and Eq. (22). As in the [integer case](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE#51-sums), it is possible to estimate these nested sums using counting techniques, for which the products of Eq. (18) can be understood as combinations of $m$ elements from a total of $n$ of the set $\Omega = \set{\lambda_1^{\alpha},\lambda_2^{\alpha},\cdots,\lambda_n^{\alpha}}$. Therefore it is necessary to use the following operations:
- [x] **nchoosek(L,m)**: generates a matrix whose rows are the different combinations of $m$ elements of the total contained in the L vector.
- [x] **prod(Array,2)**: makes the product of the arrays contained in each row of the matrix generated by the nchoosek(L,m) function.
- [x] **sum(vec)**: computes the sum of the product of the last step.

The following code includes the last instructions:

### Code #1 ###
```MATLAB
%Sum with replacement
function S1 = Suma(m,L)
COMB = nchoosek(L,m);
Vec = prod(COMB,2);
S1 = sum(Vec);
end
```
### 5.4 Shifted sums
As it was shown in Eq. (23), the sum $S_{m,n}^i$ can be computed with the same algorithm used for the sum without restriction given by $S_{m,n}$. It is only necessary to remove the term that corresponds to the $i-$ position. The following MATLAB function carries out these operations:

### Code #2 ###
```MATLAB
%Sum with replacement
function S1i = Sumai(i,m,L)
COPY = L;
COPY(i)=[];
S1i=Suma(m,COPY);
end
```
> [!IMPORTANT]
> The **Code 2** depends on the **Code 1**, and therefore it must to appear before this last one. 

### 5.5 Polynomials.
### 5.5.1 Poly_Coeff
The function **Poly_Coeff** computes the coefficients of the $P_f(s)$ polynomial given in Eq. (15). It returns a vector with real numbers and admits the following arguments:
1. A vector **L_f** that contains the lambda constants with the $\alpha$ power, **L_f**=($\lambda_1^\alpha,\lambda_2^\alpha,\cdots,\lambda_n^\alpha$)
2. A parameter **LAM_f** whose value is the same that $\Lambda^\alpha$.
3. The reactivity, denoted by the variable $\rho$.
4. A vector called **Betas** that contains the fractions of the precursors of the delayed neutrons, given by Betas =($\beta_1,\beta_2,\cdots,\beta_n$)
The following code contains the **Poly_coeff** function:
### Code #3 ###
```MATLAB
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
```
### 5.5.2 Poly_Coeff_d
The function **Poly_coeff_d** builds the coefficients of the derivative of the polynomial $P_f(s)$. It requires the vector generated by the **Poly_coeff**, and uses the following relationship:

$$ \frac{dP_f(s)}{ds}= \sum_{k=1}^{K}c_k ks^{k-1}, \tag{29}$$ 

where $c_k$ are the coefficients of the polynomial $P_f(s)$.The following code implements the Eq. (29) and returns a vector with the coefficients of $\frac{dP_f}{ds}$.
### Code #4 ###
```MATLAB
function P_d = Poly_Coeff_d(P)
exp = size(P,2)-1:-1:1;
P(1:size(P,2)-1);
P_d(1:size(P,2)-1)=P(1:size(P,2)-1).*exp;
end
```
### 5.5.3 Poly_Coeff_H
The function **Poly_Coeff_H** returns the coefficients of the polynomial $H_f(s)$ given in Eq. (16). It function requires the same input arguments that the polynomial **Poly_Coeff**, but in addition it needs a vector, denoted by C_0, whose content are the initial conditions of the precursors of the delayed neutrons. 
The following code implements the **Poly_Coeff_H**
### Code #5 ###
```MATLAB
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
```
### 5.5.4 Poly_Coeff_Q
The function **Poly_Coeff_Q** computes the coefficients of the Polynomial $Q_f(s)$ given in the Eq. (17). This polynomial only requires the vector **L_f** as input. The following code implement this function the MATLAB programming language:
### Code #6 ###
```MATLAB
function Q1 = Poly_Coeff_Q(L_f)
C_Q = [ ];
C_Q(1)=1;
for j=1:size(L_f,2)
    C_Q(j+1)=Suma(j,L_f);
end
Q1=C_Q;
end
```
## 6. AnalyticFNPKE_Insertion.m
The code **AnalyticFNPKE_Insertion.m** that is provided in the repository, solves the system given in Eq. (1) and Eq. (2), considering a constant reactivity and including the codes that were discussed before. 

>[!WARNING]
> The AnalyticFNPKE-Insertion.py code only can be used for cases with constant reactivities. For linear-time reactivities see the AnalyticFNPKE-Ramp.py code.

### 6.1 Verification of the AnalyticFNPKE_Insertion.m
A first step in the development of the codes consists of reproducing the integer case for $\alpha=1$. In other words, the fractional model must to reproduce the data obtained with the integer code [AnalyticNPKE-Insertion.py](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE#6-analyticnpke-insertionpy). Therefore we will use the following initial conditions:

$$ n (0) =1  \ \ \ \ C_i(0)=\frac{\beta_i n_0}{\lambda_i^\alpha \Lambda^\alpha}, \tag{30} $$

whose justification is provided in the submitter paper. Similarly, the following data will be used:

|Nuclear parameter | Value  ($\mathrm{s^{-1}}$)| Nuclear parameter | Value           |
| ------------- | ------------- | -------------     | --------------  |
| $\lambda_1$   |0.0127         | $\beta_1$         |0.000285         |
| $\lambda_2$   |0.0317         | $\beta_2$         |0.0015975        |
| $\lambda_3$   |0.115          | $\beta_3$         |0.00141          |
| $\lambda_4$   |0.311          | $\beta_4$         |0.0030525        |
| $\lambda_5$   |1.40           | $\beta_5$         |0.00096          |
| $\lambda_6$   |3.87           | $\beta_6$         |0.000195         |

with $\beta=0.0075$ and $\Lambda=0.0005 \mathrm{s}$. A negative reactivity given by $\rho=-1$ dollar will be used as well as a time of $t=10$ seconds. The following part of the code is considered as the "Input" section. 

### Input:

<details><summary>CLICK HERE to expand the input of the application of the AnalyticFPKE_Insertion.m</summary>
<p>

```MATLAB
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
```
</p>
</details>

> [!WARNING]  
> The input contains two different variables for the time: Target and step. This last is due to numerical issues that the Mittag-Leffler can face for very large times. Therefore, we suggest to consider that these variables have the same value for times lower than 10s. Nevertheless for greater times we suggest to use a time step, which will provide a more efficient implementation.

### Output
The output for the data provided before is give as follows:
```MATLAB
Columns 1 through 6

                        10         0.236110650788775          41.2716360082132          82.2178257032565          12.8950840417643          6.09715513377952

  Columns 7 through 8

          0.33542103625059        0.0240863276763245
```
where the first value is the Target time, the second is the neutron density and the rest of the values is the concentration of each group of the delayed neutrons.The data provided before coincides with the one reported by Nahla (2010, 1626) for the neutron density. 

### 6.2 Calculations for the fractional order $\alpha = 0.9$
As a second example we will solve the previous case, but considering a fractional order different from 1. Essentially we have the same input that the previous case with the following modification:

```MATLAB
%************************** Input **********************************

%fractional order
alpha_f = 0.9;
```
The corresponding output is the following:
### Output
```MATLAB
Columns 1 through 6

                        10         0.245804709487946          12.2037243895045          26.4909471904817          4.96488165211542          2.88018630226583

  Columns 7 through 8

         0.175643889648127         0.013646526326084
```
We can observe that the neutron density obtained with the fractional order is greater than the integer case. Nevertheless, it is possible to observe that the concentration of the precursors of the delayed neutrons is lower in each case. This can be explained by the fractional lambda constants $\lambda_i$, which have a lower value than the integer case, and therefore they are decaying in slower way. 

## 7. AnalyticFNPKE_Ramp.m
As it can be explained in the repository for the [integer case](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE#7-analyticnpke-ramppy), it is possible to use the developed solution for non-constants reactivities, dividing the time in subintervals and updating the initial conditions and the reactivity in each step. This last parameter can be estimated as follows:

$$\bar{\rho}=\frac{\ \rho\left(t_n\right)+\rho(t_{n-1})}{2}. \tag{31}$$

The code **AnalyticFNPKE_ramp.m** implements the last methodology using an iterative process and applying the **AnalyticFNPKE_insertion.m** in each step. 
> [!IMPORTANT]  
> In a strictly sense we are approximating the solution for the non-constant reactivity, which mainly depend on the size of the small step. In fact, the smaller it is the better approximation. We suggest using $h=0.0001$, which guarantees
> 
### 7.1 Verification of the AnalyticFNPKE_Ramp.m
As in the previous case, the first example of application will consists of reproducing the integer's solution assuming that $\alpha=1$. A time step of $h=0.0001 \ \ \mathrm{s}$ will be used with the purpose to have an accuracy of at least seven digits, comparing with the data reported by Hamada (2018, p. 7). The following table contains the values of the decay constants as well as the fractions betas proposed by Nahla (2010, p. 1626):

|Nuclear parameter | Value  ($\mathrm{s^{-1}}$)| Nuclear parameter | Value           |
| ------------- | ------------- | -------------     | --------------  |
| $\lambda_1$   |0.0127         | $\beta_1$         |0.000266         |
| $\lambda_2$   |0.0317         | $\beta_2$         |0.001491         |
| $\lambda_3$   |0.115          | $\beta_3$         |0.001316         |
| $\lambda_4$   |0.311          | $\beta_4$         |0.002849         |
| $\lambda_5$   |1.40           | $\beta_5$         |0.000896         |
| $\lambda_6$   |3.87           | $\beta_6$         |0.000182         |


The initial conditions are the same that in the previous case, $\beta=0.007$ and $\Lambda= 0.00002 \ \mathrm{s}$ and the reactivity is given by:
$$\rho(t)=0.1\beta t, \tag{32}$$

### Input
The following lines show how the input parameters can be introduced to the code.

<details><summary>CLICK HERE to expand the input of the application of the AnalyticFPKE_ramp.m</summary>
<p>

```MATLAB
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
```
</p>
</details>

### Output ###
The **AnalyticFNPKE_ramp.m**  generates a single output, which consits of ".xlsx" file with a single sheet in which the first column is the time, the second one is the neutron density and the rest of columns contains the results of the precursors of delayed neutrons. 

### Output ###
The output file for this example is given in the following [link](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticalFNPKE/blob/main/Output_ramp.xlsx).

> [!WARNING]  
> The size of the .xlsx is about 10-20 Mb.

> [!NOTE]  
> The code provided a counter where is showed the progress of the calculation, which appears in the command Window of the MATLAB software.

### Verification ###
The following table contains a comparison of the data computed (extracted from the .xlsx) with the data reported by Hamada (2018, p. 3034). As it can be observed, the results agree with at least seven precision digits. 

|Time (s) | Reference (Hamada, 2018)| AnalyticFNPKE_ramp | Relative %          |
| ------------- | ------------- | -------------     | --------------  |
| 2.0           |1.3382000      | 1.3382000         |1.10394E-08      |
| 4.0           |2.2284418      | 2.2284418         |3.18368E-08      |
| 6.0           |5.5820524      | 5.5820522         |1.14074E-07      |
| 8.0           |42.786295      | 42.786294         |5.16193E-07      |
| 10.0          |451163.62      | 451163.61         |0.009243646      |

As it can be observed the error remains below of $10^{-3}$ %. 

### 7.2 Calculations for the fractional order $\alpha = 0.9$
As a second example, the **AnalyticFNPKE_ramp.m** code will be used to compute the neutron's density considering a fractional order $\alpha =0.5$. The input is very similar to the one used in the past case, but the following line must to be changed:

```MATLAB

%*******************************************************************
%************************** Input **********************************
%fractional order 
alpha_f = 0.9;
```
