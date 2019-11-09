clc
clear all
close all

rho=996;%Density of water kg m^-3
Re=10000;%Re range for vortex street generation
v=0.798e-3/rho;%Kinematic Viscosity of Water @ 30 degC
D=2e-2;%Circular Cylinder Diameter/m
St=0.22;%St for circular cylinder
 
V=(Re*v)/D;%Flow speed in test section/m s^-1
fv=(2*pi)^-1*(St.*V)/D;%Vortex shedding frequency/Hz

L_f=60e-3;%Length of Fluid in Contact
t=0:0.001:1;

%C_L=0.045+1.05*(1-Re/1600)^4.5;%For 260<Re<1600
%C_L=0.045+3*(log10(Re/1600))^4.6;%For 1600<Re<5400
C_L=0.52-0.06*(log10(Re/1600))^-2.6;%For 5400<Re<225000
A=0.5*C_L*L_f*rho*D*V^2
f_v=St*V/D
w_v=2*pi*f_v
T_v=1/f_v;
Lift=A.*cos(w_v.*t);

plot(t, Lift)
title('Lift force variation with time')
xlabel('Time / (s)')
ylabel('Lift Force / (N)')
