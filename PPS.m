clc
clear all 
close all
format long

% Specific Enthalpies at respective states [kJ/kg]
h = 

h1 = 3450.47;
h2 = 2065.31;
h3 = 191.812;
h4 = 206.9;
h5 = 727.143;

vf = 0.00165696; %Specific Volume of Sat. Liquid @ 150 bar
vg = 0.0103401; %Specific Volume of Sat. Vapour @ 150 bar
hf = 1610.15; %Specific Enthalpy of Sat. Liquid @ 150 bar
hg = 2610.86;%Specific Enthalpy of Sat. Vapour @ 150 bar

hfg = hg - hf;%Saturation Length on P-h curve for 150 bar

wT = h1 - h2;%Turbine Work Out [kJ/kg]
wP = h4 - h3;%Pump Work In [kJ/kg]

W = 250e3;%Net Power Output [Given Parameter] [kW]

m = W/(wT-wP)%Mass Flow Rate of steam [kg/s]

C = 25e3;%Calorific Value of Coal [kJ/kg]
eta_stg = 0.9;%Steam Generator Efficiency

mf = m*(h1-h5)/(eta_stg*C)%Rate at which Coal is Burnt [kg/s]

ef = m/mf;%Evaporation factor

H = 45;%Height of boiler [m]
g = 9.81;%Acceleration due to gravity [m/s^2]

rhoD = 1/vf;%Downcomer (Sat. Liquid) Fluid Density

CR = [6:26]';%Range of Acceptable Circulation Ratios
xtop = 1./CR;%Corresponding possible Top Dryness Fractions

figure(5)
plot(CR,xtop)
title('TDF vs CR')
xlabel('Circulation Ratio')
ylabel('Top Dryness Fraction')
grid on
grid minor

vfg = vg - vf;%Saturation Length on P-v curve for 150 bar

vtop = vf+xtop.*vfg;%Specific Volume at the top of the Riser tubes
rho_top = 1./vtop;%Correspondig range of possible density at the top of riser tubes

rhom = 0.5.*(rhoD + rho_top);%Correspondig range of possible mixture density in riser tubes

p = g*H.*(rhoD - rhom);%Pressure Head availabe for Natural Circulation

S = 1.2; %Slip ratio

psi = S*vf/vg%Psi constant
a = 1./(1+((1-xtop)./xtop).*psi);%Void Fraction at the top of riser tubes

figure(6)
plot(CR,a)
title('\alpha vs CR')
xlabel('Circulation Ratio')
ylabel('Void Fraction')
grid on
grid minor


Dr = [62.5:0.1:76.5]'./1000;%Range of Riser (Outer) pipe dimaeters

CR = [6:26]';%Range of Circulation Ratios

figure(4)
plot(CR,p)
title('Pressure Head vs CR')
xlabel('Circulation Ratio')
ylabel('Pressure Head [Pa]')
grid on
grid minor

i = 1;

while i <= length(Dr)
    CRn(i,:) = CR;
    i = i+1;
end 

CR = CRn;
%Dr = Drn;

xtop = 1./CR;%Possible top dryness fractions


vfg = vg - vf;%Redefining Saturation Length on P-v curve for 150 bar for Dr
vtop = vf+xtop.*vfg;%Redefining Specific Volume at the top of the Riser tubes
rho_top = 1./vtop;%Redefining Correspondig range of possible density at the top of riser tubes


tr = 3e-3;%Tube wall thickness [m]
Uc = 1.8;%Circulation Velocity [m/s]

k = (rho_top.*xtop);%Term to be used for Parametric Sweep
i = 1;

while i<= 21
  
    msfr(:,i) = 0.25*pi.*(Dr-2*tr).^2.*k(:,i);%Steam formation rate at riser exit
    i = i+1;
end

nr = m./msfr; %Number of required riser tubes

i = 1;

figure(1)
while i<= 21
  plot(Dr,nr(:,i))
  hold on
  i = i+1;
end
title('Number of tubes vs Riser Diameter for CR = 6-26')
xlabel('Riser Tube Diameter [m]')
ylabel('Number of Tubes')
grid on
grid minor

i = 1;

while i<= 21
  
    Qar(:,i) = msfr(:,i).*hfg./(H.*Dr);%Steam formation rate at riser exit
    i = i+1;
end

i = 1;
figure(2)
while i<= 21
  plot(Dr,Qar(:,i))
  hold on
  i = i+1;
end
title('Heat Rate vs Riser Diameter for CR = 6-26')
xlabel('Riser Tube Diameter [m]')
ylabel('Heat Rate [kW/m^2]')
grid on
grid minor

i = 1;
figure(3)
while i<= 21
  plot(Dr,msfr(:,i))
  hold on
  i = i+1;
end
title('Rate of Steam Exiting One Riser vs Riser Diameter for CR = 6-26')
xlabel('Riser Tube Diameter [m]')
ylabel('Steam Flow Rate [kg/s]')
grid on
grid minor

Ddc = [150:250]./1000;%Range of Downcomer pipe dimaeters
mD = 0.25*pi.*(Ddc).^2.*rhoD.*Uc;%Saturated water mass flow rate in a downcomer tube

figure(7)
plot(Ddc,mD)
title('Mass Flow Rate in One Downcomer Tube vs Downcomer Tube Diameter')
xlabel('Downcomer Tube Diameter [m]')
ylabel('Mass Flow Rate [kg/s]')
grid on
grid minor

nD = m./mD; %Number of downcomer tubes

figure(8)
plot(Ddc,nD)
title('Number of tubes vs Downcomer Tube Diameter')
xlabel('Downcomer Tube Diameter [m]')
ylabel('Number of Tubes')
grid on
grid minor

hffw = 719.206;%Specific Enthalpy of Entering Feed Water (approx. Sat. Liquid @170oC) 
Tfwi = 170;
Tfwe = 342.158;

Qeco = m*(hf - hffw);%Heat Transfer in the Economiser section

mfl = 1800;
cp = 1.12;%Specific Heat at Const. Pressure of flue gas [kJ/kg.K]

Tfe = 450;% Flue Gas Exit Temperature in Economiser Section [oC]

Tfi = Qeco/(mfl*cp)+Tfe% Inlet Flue Gas Temperature [oC]

dTi = Tfi - Tfwi;
dTe = Tfe - Tfwe;
Tlm = (dTi - dTe)/log(dTi/dTe)

Do = 70e-3;%Economiser Tube Outer Diameter
te = 5e-3;%Economiser Tube thickness
Di = Do - 2*te;%Economiser Tube Inner Diameter

Uo = 80e-3;%Overall Heat Transfer Coefficient [kW/m^2.K]

Ao = Qeco/(Uo*Tlm);%conomiser Total Heat Exchange Surface

meco = 0.25*pi*Di.^2*Uc/vf;

neco = m./meco;
neco = 60;

leco = Ao/(neco*pi.*Do);

B = 5;
C = 5e-3;
nt = leco./(B-2*C);
nt = 148;

pitch = 45e-3;

Heco = pitch.*nt;

ma = 1500;
Tai = 35;

cp = 1.12;%Specific Heat at Const. Pressure of flue gas [kJ/kg.K]
cpa = 1.001;%Specific Heat at Const. Pressure of air [kJ/kg.K]

Tfe = 150;% Flue Gas Exit Temperature in Economiser Section [oC]

Tfi = 450;% Inlet Flue Gas Temperature [oC]

Qph = mfl*cp*(Tfi-Tfe);%Heat Transfer in the Air Pre Heater section

Tae = Qph/(ma*cpa)+Tai;

dTi = Tfi - Tae;
dTe = Tfe - Tai;
Tlm = (dTi - dTe)/log(dTi/dTe)

Do = 65e-3;%PH Tube Outer Diameter
Di = 60e-3;%PH Tube Inner Diameter

Uo = 30e-3;%Overall Heat Transfer Coefficient [kW/m^2.K]

Ao = Qph/(Uo*Tlm);%Air PH Total Heat Exchange Surface

R = 0.287;
vfl = R*(Tfi+273)/101.325;
Uf = 12;

mph = 0.25*pi*Di.^2*Uf/vfl;

nph = mf/mph
nph = 1334;

lph = Ao/(nph*pi.*Do);

hi = 2610.86;

Qsh = m*(h1-hi)
qm = 140;
Do = 60e-3;
t = 5e-3;
Di = Do - 2*t;
Ush = 10;
msh = 0.25*pi*Di.^2*Ush/vg;

nsh = m/msh;
nsh = 97;

Ao = Qsh/qm;

lsh = Ao/(nsh*pi.*Do)









