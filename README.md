# Analytical Solution of the Fractional Neutron Point Kinetic Equations. 
The present repository contains the MATLAB codes developed to solve the Fractional Neutron Point Kinetic Equations (FNPKE). These codes were described in the article *A New Solution of the Fractional Neutron Point Kinetics Equations using Symmetry and the Heaviside’s expansion formula*, which was recently submitted to the **Progress in Nuclear Energy** journal. 

The programs are licensed under a Creative Commons Attribution 4.0 International License: http://creativecommons.org/licenses/by/4.0/

Authors: Carlos-Antonio Cruz-López (cacl.nucl@gmail.com), Gilberto Espinosa-Paredes (gepe@xanum.uam.mx)

Mathematical and algorithmical generalities of the codes are described in the following lines with the purpose to provide some insight of the developed work. Nevertheless, a more detailed discussion is provided in the submitted article.
## Software specifications.
The AnalyticFNPKE codes were written in the MATLAB programming language R2021a, but older versions can be used, from the R2012b one onward. The reported examples and results were obtained in a 3.8 GHz desktop computer, with 32 Gb of RAM and under a Windows 11 environment. The codes require the following script written by Roberto Garrappa (2015, 2024):
- [x] ml.m

This last file can be freely download from the mathworks site in the following [link](https://www.mathworks.com/matlabcentral/fileexchange/48154-the-mittag-leffler-function).

## Financial Support.
The authors appreciate the financial support received from the Consejo Nacional de Humanidades, Ciencia y Tecnología (CONAHCYT), under the program *Estancias Posdoctorales por México, 2022*, with the project entitled: *Desarrollo de modelos fenomenológicos energéticos de orden fraccional, para la optimización y simulación en reactores nucleares de potencia*, by which the present development was possible.

## 1. Mathematical description of the problem
Following the ideas of Ray and Patra (2014), a compartmental version of the The Fractional Neutron Point Kinetic Equations (FNPKE), with $K$ groups of precursors of delayed neutrons, can be written as follows:

$$D_C^\alpha n\left(t\right)=\frac{\rho\left(t\right)-\beta}{\Lambda^\alpha}n\left(t\right)+\sum_{i=1}^{K}{\lambda_i^\alpha C_i(t)}, \tag{1}$$
