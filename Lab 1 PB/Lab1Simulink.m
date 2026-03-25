% clear
% close all
% clc

Tau1 = 40*10^-3; % ms => s
Tau2 = 40*10^-3;
J = 0.00276; % N*ms^2/rad
G = 0.03; % N*ms/rad
K = 0.992; % N*m/rad

a4 = Tau1*Tau2*J;
a3 = (Tau1*Tau2*G+J*Tau1+J*Tau2)/a4;
a2 = (Tau1*Tau2*K+G*Tau1+G*Tau2+J)/a4;
a1 = (K*Tau1+K*Tau2+G)/a4;
a0 = K/a4;

b0 = 1/a4;

