clc
clear 
format long 

xc=[0:0.125:0.5]';
xh=[0 0.5]';

%C1
Th1=[54 49]';
Tc1=[39 44 44 44 45]';

%C2
Th2=[55 49]';
Tc2=[47 46 44 43 40]';

plot(xh, Th2, xh, Th1, xc, Tc2, xc, Tc1) 

legend('Configuration 2: Hot', 'Configuration 1: Hot', 'Configuration 2: Cold', 'Configuration 1: Cold')

%C3
Th3=[51 47]';
Tc3=[44 43 43 41 39]';

%C4
Th4=[48 45]';
Tc4=[43 42 41 40 38]';

%C5
Th5=[44 42]';
Tc5=[41 41 40 39 37]';

Th=[Th2 Th3 Th4 Th5]'
Tbh=(Th(:,1)+Th(:,2))./2%Bulk Temperatures for Hot Stream
Tc=[Tc2 Tc3 Tc4 Tc5]'
Tbc=(Tc(:,1)+Tc(:,5))./2%Bulk Temperatures for Cold Stream

Vh=[2 2 2 2]';
Vc=[2 2 2 2]';

Vh=0.001*60^-1*Vh%Volumetric flow rate in m^3 s^-1
Vc=0.001*60^-1*Vc%h:hot c:cold

%Ah=21*0.25*pi*0.008^2%Effective cross sectional area for hot water flow
%Ac=(0.25*pi*0.148^2)-(21*0.25*pi*0.01^2)%Cross sectional area for cold water flow

rhoh=[986.94 988.5 989.5 990.9]';%Densities for hot water evaluated at or near bulk temperatures
rhoc=[990.7 991.5 991.9 992.48]';%Densities for cold water evaluated at or near bulk temperatures

mh=rhoh.*Vh%Mass flow rate of hot water 
mc=rhoc.*Vc%Mass flow rateof cold water

cph=[4181.8 4180.8 4180.3 4179.6]';%Specific Heats for hot water evaluated near bulk temperatures
cpc=[4176.7 4179.3 4179.1 4178.8]';%Specific Heats for cold water evaluated near bulk temperatures

qh=mh.*cph.*(Th(:,1)-Th(:,2))
qc=mc.*cpc.*(Tc(:,1)-Tc(:,5))

ql=qh-qc

delT1=Th(:,1)-Tc(:,1)
delT2=Th(:,2)-Tc(:,5)

delT_lm=(delT1-delT2)./(log(delT1./delT2))%Logarithmic Mean Temperature difference

U=qh./(0.00192.*delT_lm)

Ttank=[60 55 50 45]';

plot(Ttank, U);






