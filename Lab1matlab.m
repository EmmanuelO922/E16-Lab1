% clear
% close all
clc

%% Analytical Solution

syms theta(t) ui(t)

Tau1 = 40*10^-3; % ms => s
Tau2 = 40*10^-3;
J = 0.00276; % N*ms^2/rad
G = 0.03; % N*ms/rad
K = 0.992; % N*m/rad

ui(t) = heaviside(t) - heaviside(t - 1); 

Dtheta = diff(theta,t);
D2theta = diff(theta,t,2);
D3theta = diff(theta,t,3);
F = J*diff(theta,t,2)+G*diff(theta,t)+K*theta; % + 0.1;
% I just defined what F was in terms of theta and instead of solving for f
% I plugged it into the ODE and solved for theta
ode = Tau1*Tau2*diff(F,t,2)+(Tau1 + Tau2)*diff(F,t)+F == ui;
conds= [theta(0)==0, Dtheta(0)==0, D2theta(0)==0, D3theta(0)==0];
thetaSol = dsolve(ode, conds);
disp('Solution to the ODE = ')
disp(thetaSol)

figure; hold on
fplot(ui,[0,2],'k',LineWidth=0.75)
fplot(thetaSol,[0,2],'b',LineWidth=1.5)

%% Numerical Solution


[t,y] = ode45(@jointang,[0,2],[0;0;0;0]);
plot(t,y(:,1),'--','Color','r',LineWidth=1.5)


%% Functions

function dx = jointang(t,x)

Tau1 = 40*10^-3;
Tau2 = 40*10^-3;
J = 0.00276;
G = 0.03;
K = 0.992;

a4 = Tau1*Tau2*J;
a3 = (Tau1*Tau2*G+J*Tau1+J*Tau2)/a4;
a2 = (Tau1*Tau2*K+G*Tau1+G*Tau2+J)/a4;
a1 = (K*Tau1+K*Tau2+G)/a4;
a0 = K/a4;
b0 = 1/a4;

f = heaviside(t) - heaviside(t - 1);
dx1 = x(2);
dx2 = x(3);
dx3 = x(4);
dx4 = -a3*x(4)-a2*x(3)-a1*x(2)-a0*x(1)+b0*f;

dx = [dx1 dx2 dx3 dx4]';

end