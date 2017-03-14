clear all;
clf

P0 = 343233 %Pressure
Rg = 8.314  %gas constant
T0 = 823  %Initial temperature
Ea = [48906]  %Activation Energy
A0 = 3715.35 %Pre Exponential Factor
k0 = A0*exp(-Ea/(Rg*T0)) %Rate Constant
R =  1.85/2  %outer radius of the reactor
r = 0.551  %inner radius of the reactor
K = r/R
VR0 = 3.34*10^(-1) %maximum radial velocity
Mw = 94/1000  %Molecular weight of the component
Ci0 = P0/(Rg*T0)  %maximum concentration of the component
rho0 = 4.92 %density of mixture
H0 =  204870  %molar enthalpy of reaction
Cp =  Rg*(3.148+18.4*(10^(-3))*T0+1.36*(10^(-4))*(T0^2)-1.88*(10^(-7))*(T0^3)+7.36*(10^(-11))*(T0^4)) %molar heat capacity

S = [-1 1 3] % p X N
Sb = [-1] % p x 1
chi = abs(sign(S)) % p X N
SumS = sum(S,2)
F0 = [784225 0 0]

A = (k0*R)/VR0
C = Ea/(Rg*T0)
D = H0/(Cp*T0)
G = R*VR0/P0
I = R*VR0*VR0/P0
J = rho0*VR0/R
Er1 = 5000
Er2 = 30000

tspan = [1 0];

%    phi_A phi_B phi_C theta P Vr rho0
y0 = [1      0     0     0   1 -1  4.92];
yp0 = [0; 0; 0; 0; 0; 0];
[y0_new,yp0_new] = decic(@model,tspan,y0,[1 1 1 1 1 1],yp0,[]);

myrxnfun = @twofirst;
N = 4;

f = @(t,y,yp) model(t,y,yp,myrxnfun,N,Er, G, I , J)

[t,y] = ode15i(f,tspan,y0_new,yp0_new);



plot(t,y(:,1))
title('phi_A')
figure
plot(t,y(:,2))
title('phi_B')

figure
plot(t,y(:,3))
title('phi_C')

figure
plot(t,y(:,4))
title('Temperature')
figure

plot(t,y(:,5))
title('Pressure')

figure
plot(t,y(:,6))
title('Radial velocity')
