clc
clear 

L= 0.1;%Plate Length (m)

w=0.07; %Plate width (m)

pt=0.015;%Plate thickness (m)

A=2*((L*w)+(L*pt)+(w*pt)); %Surface Area (m^2)

V=L*w*pt; %Volume (m^3)

Lc=V/A %Characteristic length(m)

t=[0:10:70]';

a=0.45e-3; %Thermal diffusivity (m^2 s^-1)

Fo=a.*t./Lc^2

T1=[73.1 73.4 73.8 74 74.3 74.7 75 75.3]';

T2=[28.9 56.7 59.6 61.1 62.3 63.7 65 66.1]';

T3=[27.9 29.6 37.3 44.9 49.9 54.7 57.9 60.8]';

rho=8500; %Density (kg/m^3)
k=16.3 %Thermal conductivity (W/m.K)
c=460 %Specific heat capacity (J/kg)

Bi=0.1

h=Bi*k/Lc%Convection heat transfer coefficient

Tau=(rho*Lc*c)./h

T3_t=exp(-t./Tau).*(T3(1)-T1)+T1;

plot(t, T3_t, t, T3)
grid on
legend('Theoretical', 'Experimental')

Qmax=rho*V*c.*(80-T3(1));

Q=h*A.*(T2-T1)

plot(t, Q)
grid on
