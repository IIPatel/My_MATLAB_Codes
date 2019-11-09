clc
clear
format long

Re=80:200;%Re range for vortex street generation
v=0.798e-3/996;%Kinematic Viscosity of Water @ 30 degC
D=1e-2;%Circular Cylinder Diameter/m
St=0.22;%St for circular cylinder
fn=2.25;%Natural Frequency of cylinder/Hz


%Vr=V./(fn.*D);%Reduced velocity
%fv=(St.*V)./D;%Vortex shedding frequency/Hz

A=(25e-2)*(20e-2);% Test section flow CS Area/m^2
%Re=(Dh.*V)./v
Dh=(4*A)/(2*(20e-2)+(25e-2));
V=(Re.*v)./Dh;%Require flow speed in test section/m s^-1

Q=A.*V;%Required volume flow rate/m^3 s^-1
Q=Q.*(1000*60);%L/min
subplot(221)
plot(Re, Q)
xlabel('Reynold Number')
ylabel('Volumetric Flow rate/ L min^-1')
% subplot(222)
% plot(Re,fv)
% xlabel('Reynold Number')
% ylabel('Vortex Shedding Frequency / Hz')

% o=250e-3/2;%Outlet half-width
% i=o*2;%Inlet half-width;contraction ratio is i/o
% L=600e-3;%Length of contraction section
% x=[0:1e-3:L]';%x coordinate points of section
% y=i-(i-o).*(6.*(x./L).^5-15.*(x./L).^4+10.*(x./L).^3);%y coorinates to model 5th order curve
% subplot(223)
% plot(x,y)

% filename = 'Contraction Section.xlsx';
% sheet = 1;
% xlRange = 'A2';
% xlswrite(filename,x,sheet,xlRange)
% filename = 'Contraction Section.xlsx';
% sheet = 1;
% xlRange = 'B2';
% xlswrite(filename,y,sheet,xlRange)