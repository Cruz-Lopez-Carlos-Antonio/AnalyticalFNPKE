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

$$\lim_{\alpha \rightarrow 1} \underbrace{S_m(\Omega)}_ {\mathrm{Fractional\ case}}= S_m(\lim_{\alpha \rightarrow 1}{\Omega})=S(\lim_{\alpha \rightarrow 1} \set{\lambda_1^\alpha,\lambda_2^\alpha,\cdots,\lambda_n^\alpha}) = S_m(\set{\lambda_1,\lambda_2,\cdots,\lambda_n})

It is possible to write the Polynomials as two variable functions that depend on $s$ as well as on $\alpha$:


=S\left(\lim_{\alpha \rightarrow 1}{\Omega}\right)=S\left(\lim_{\alpha\rightarrow1}{\left\{\lambda_1^\alpha,\lambda_2^\alpha,\cdots,\lambda_n^\alpha\right\}}\right)


