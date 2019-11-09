clc
clear all 
close all

rho_bf=2864.7889756541160438399077407053;%Bluff Body material density / kg m^-3
b=1e-2;%Major axis half length in m
Lbf=6e-2;%Length of Bluff Body in m

k=68e-3;%Radius of Gyration=Lbi

r=[1 4 6];%Aspect ratios
a=b.*r;%Major axis half length in m
Mbf=rho_bf*pi*b^2*Lbf.*r%Mass of bluff body in kg

%Mbf.*(3.*(a./r).^2+4*Lbf^2)/12
Jbf =Mbf.*(3.*b^2+4*Lbf^2)/12;%Mass moment of Inertia of the Bluff Body about its CG ;
Jtip=Jbf+Mbf.*(0.5*Lbf)^2;%Mass moment of Inertia of the Bluff Body about the bimorph tip ;
m=Jtip/(k^2)%Pseudo Mass with CG at radius of gyration

L=;%Length of Pseudo Mass
W=1;%Width fixed
H=1.7e-3;%Height
V=L*W*H;%Volume
rho=m/V%Required Density

rho=1000;%Density of water kg m^-3
v=0.798e-3/rho;%Kinematic Viscosity of Water @ 30 degC
D=2e-2;%Circular Cylinder Diameter/m
St=0.22;%St for circular cylinder 

V=0.31;%Flow speed in test section/m s^-1
Re=V*D/v%Re #
L_f=55e-3;%Length of Fluid in Contact
t=0:0.001:1;

%C_L=0.045+1.05*(1-Re/1600)^4.5;%For 260<Re<1600
%C_L=0.045+3*(log10(Re/1600))^4.6;%For 1600<Re<5400
C_L=0.52-0.06*(log10(Re/1600))^-2.6;%For 5400<Re<225000
A=0.5*C_L*L_f*rho*D*V^2
f_v=St*V/D
w_v=2*pi*f_v;
T_v=1/f_v;
Lift=A.*cos(w_v.*t);

