clc
clear
format long
 
Re=100:10:10000;%Re range for vortex street generation
v=0.798e-3/996%Kinematic Viscosity of Water @ 30 degC
D=2e-2;%Circular Cylinder Diameter/m
St=0.22;%St for circular cylinder
 
V=(Re.*v)./D;%Flow speed in test section/m s^-1
fv=(St.*V)./D;%Vortex shedding frequency/Hz
 
A=(25e-2)*(25e-2);% Test section flow CS Area/m^2
 
Dht=(4*A)/(2*(20e-2)+(25e-2))%Hydraulic Diameter of test channel
 
Q=A.*V%Volumetric flow rate/m^3 s^-1
Q=Q.*(1000*60)%Conversion into L/min
subplot(221)
plot(Re, Q)
title('Q against Re')
xlabel('Reynolds Number')
ylabel('Volumetric Flow rate / (L/min)')
 
subplot(222)
plot(Re,fv)
title('f_v against Re')
xlabel('Reynolds Number')
ylabel('Vortex Shedding Frequency / (Hz)')
 
o=250e-3/2;%Outlet half-width
i=o*2;%Inlet half-width; contraction ratio is i/o
Lc=600e-3;%Length of contraction section
xc=[0:1e-3:Lc]';%x coordinate points of section
yc3=i-(i-o).*(-2.*(xc./Lc).^3+3.*(xc./Lc).^2);%y coordinates to model 3rd order curve
yc5=i-(i-o).*(6.*(xc./Lc).^5-15.*(xc./Lc).^4+10.*(xc./Lc).^3);%y coordinates to model 5th order curve
yc7=i-(i-o).*(-20*(xc./Lc).^7+70.*(xc./Lc).^6-84.*(xc./Lc).^5+35.*(xc./Lc).^4);%y coordinates to model 7th order curve
zc=zeros(601,1);%z coordinates 
subplot(223)
plot(xc,yc5)
title('Contraction Section (2:1)')
ylabel('z / (m)')
xlabel('x / (m)')
hold on
plot(xc,yc3)
hold on
plot(xc,yc7)
legend('3rd order contour', '5th order contour','7th order contour')
 
Lt=0.8;%Length of contraction section
xt=[0.6:1e-3:(Lc+Lt)]';%x coordinate points of section
for i=1:length(xt)
    yt(i)=o;%y coordinates
    i=i+1e-3;
end
yt=yt';
zt=zeros(801,1);%z coordinates to model 5th order curve
plot(xt,yt)
subplot(224)
plot(Re,V)
title('V against Re')
xlabel('Reynolds Number')
ylabel('Flow Velocity / (m/s)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lbi = 76.68e-3; %Length of the overhung piezobeam in m
Lp = Lbi-1.5e-2; %Approximate Length of PZT layers
W = 32.98e-3-0.05e-3;%Width of the piezobeam in m
Hs = 0.22e-3;%Thickness of the shim in m
Hp = 0.22e-3; %Thickness of the piezo in m
Es = 11.2e10; sm = 1/Es; %Young’s Modulus and compliant of the shim in Pa and Pa^-1
Ep = 5.6e10; s11 = 1/Ep; %Young’s Modulus and compliant of piezo in Pa and Pa^-1
rho_s = 8780; %Density of shim in kg m^-3
rho_p = 7500; %Density of piezo in kg m^-3

mass_beampiezo = 2*rho_p*Lp*W*Hp;
mass_beamshim = rho_s*Lbi*W*Hs;
Mbi = mass_beampiezo + mass_beamshim; %Total mass of bimorph

Ip=2*W*((Hp^3)/12+Hp*(0.5*Hp+0.5*Hs)^2);%Second moment of area about the neutral axis of the PZT layers.
Is=(W*Hs^3)/12;%Second moment of area about the neutral axis of the shim.
EI=Ep*Ip+Es*Is;%Effective Flexural Stiffness of the bimorph [Pa.m^4]

k_eq=3*EI/Lbi^3;%Equivalent Spring Stiffness of the solid structure [N/m]

%Bluff body parameters
rho_bf = [1300 1070];%Bluff Body material density [PLA ABS] in [kg/m^3]
a=1e-2;%Major axis half length in [m]
Lbf=[8e-2 10e-2];%Length of Bluff Bodies [C E] in [m]

r = [1 2];%Aspect ratios

i = 1;

while i<=length(Lbf)
    
b(i) = a./r(i)%Minor axis half length in [m]

Mbf(i) = rho_bf(i).*pi*a^2.*Lbf(i)/r(i);%Mass of bluff body in [kg]
Jbf(i) = Mbf(i)*(3.*(a/r(i)).^2 + 4*Lbf(i)^2)/12;%Mass moment of Inertia of the Bluff Bodies about the its CG [kg.m^2]
Jtip(i) = Jbf(i)+Mbf(i)*(0.5.*Lbf(i))^2;%Mass moment of Inertia of the Bluff Bodies about the bimorph tip [kg.m^2]

Lfc = [4e-2 6e-2];%Approximate Fluid Contact Length [m]  
rho_f = 996;%Fluid Density [kg/m^3]
M_f(i) = Lfc(i)*rho_f*pi.*a.^2.*b(i);%Fluid Added Mass [kg]

m_eq(i) = (33/140)*Mbi + Mbf(i) + 4.*Jbf(i)/(Lbf(i))^2 + M_f(i);%Equivalent Mass [kg]
f_s(i) = 1/2/pi.*sqrt(k_eq./m_eq(i));%Natural Un-damped Frequency [Hz]

i = i+1;

end
