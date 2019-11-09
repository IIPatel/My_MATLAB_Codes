clc
clear all
close all

Lbi= 76.68e-3; %Length of the overhung piezobeam in m
Lp=Lbi-1.5e-2;
W= 32.98e-3-0.05e-3;%Width of the piezobeam in m
Hs= 0.22e-3;%Thickness of the shim in m
Hp = 0.22e-3; %Thickness of the piezo in m
Es = 11.2e10; sm = 1/Es; %Young’s Modulus and compliant of the shim in Pa and Pa^-1
Ep = 5.6e10; s11 = 1/Ep; %Young’s Modulus and compliant of piezo in Pa and Pa^-1
rho_s = 8780; %Density of shim in kg m^-3
rho_p = 7500; %Density of piezo in kg m^-3

mass_beampiezo = 2*rho_p*Lp*W*Hp;
mass_beamshim = rho_s*Lbi*W*Hs;
Mbi = mass_beampiezo + mass_beamshim %Total mass of bimorph
mass_dist = Mbi/Lbi; % Distributed mass per length

Ip=2*W*((Hp^3)/12+Hp*(0.5*Hp+0.5*Hs)^2);%Second moment of area about the neutral axis of the piezo layers.
Is=(W*Hs^3)/12;%Second moment of area about the neutral axis of the shim.
EI=Ep*Ip+Es*Is;%Effective Flexural Stiffness of the bimorph

k_eq=3*EI/Lbi^3;%Equivalent Spring Stiffness of the solid structure

%Bluff body parameters
rho_bf=1300;%1070;%7850;%Bluff Body material density / kg m^-3
a=1e-2;%Major axis half length in m
Lbf=10e-2;%Length of Bluff Body in m
r=[1 2];%Aspect ratios
b=a./r.*100%Minor axis half length in cm
Mbf=rho_bf*pi*a^2*Lbf./r%Mass of bluff body in kg
Jbf=Mbf.*(3.*(a./r).^2+4*Lbf^2)/12;%Mass moment of Inertia of the Bluff Body about the its CG
Jtip=Jbf+Mbf.*(0.5*Lbf)^2%Mass moment of Inertia of the Bluff Body about the bimorph tip ;
rho_f = 1000;%Fluid Density
M_f = Lbf*rho_f*pi.*a.^2.*b;%Fluid Added Mass


m_eq=(33/140)*Mbi+Mbf+4.*Jbf/(Lbf)^2+M_f;%Equivalent Mass

f_s=(2*pi)^-1*sqrt(k_eq./m_eq)%Structural (Natural) Frequency / [Hz]

% k=Lbi;%Radius of Gyration / [m]
% m_ps=Jtip./(k^2);%Pseudo Mass with CG at radius of gyration
% 
% % L=4e-3;%Length of Pseudo Mass
% % H=1.7e-3;%Height
% % V=(L*H+(0.25e-3)^2);%Volume
% % rho=m_ps/V%Required Density





