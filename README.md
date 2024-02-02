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
