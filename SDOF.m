clc
clear all
close all
%Piezoelectric Coefficient

d_31= -186e-12;
d_33=6.70e-010;

%Elastic Compliance
sE_33=2.07e-011;
sE_11=1.65e-011;
sE_12=-4.78e-012;
sE_13=-8.45e-012;
sE_31=-8.45e-012;

cE_11=sE_11/(sE_11^2-sE_12^2);
cE_33=sE_33^-1; 

%Planar Coupling Coefficient
k_p=0.65;%From vendor provided material data
sigmaP = -sE_12/sE_11;%Planar Poisson Ratio
epsT_33=8.854187817e-12*3400;%Permittivity at constant stress

epsS_33 = epsT_33 - 2*d_31^2/(sE_11+sE_12);%Permitivity at constant strain
%epsS_33=epsS_33/8.854187817e-12;
e_31=d_31/(sE_11+sE_12);

%e_33=d_31/sE_11;
%epsS_33 = epsT_33 - d_31*e_33;%Permitivity at constant strain

k_31=sqrt(k_p^2*0.5*(1-sigmaP))%Coupling Coefficient for {3-1} Actuation Mode
k_e=sqrt(k_31^2./(1-k_31^2))%Effective Coupling for electromechanical coupled model
%k_e=sqrt((e_31^2./(cE_33*epsS_33))/(1-(e_31^2./(cE_33*epsS_33))))

%e_31=sqrt(k_e^2*sE_33^-1*epsS_33)

zeta_m=0.2360;%Effective Damping Ratio; Obtained from Free Vibration Tests in Fluid Medium

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lbi = 76.68e-3; %Length of the overhung piezobeam in m
Lp = Lbi-1.5e-2;
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

Ip=2*W*((Hp^3)/12+Hp*(0.5*Hp+0.5*Hs)^2);%Second moment of area about the neutral axis of the piezo layers.
Is=(W*Hs^3)/12;%Second moment of area about the neutral axis of the shim.
EI=Ep*Ip+Es*Is;%Effective Flexural Stiffness of the bimorph

k_eq=3*EI/Lbi^3;%Equivalent Spring Stiffness of the solid structure

%Bluff body parameters
rho_bf=[1300 1070];%Bluff Body material density [PLA ABS] in [kg/m^3]
a=1e-2;%Major axis half length in [m]
Lbf=[8e-2 10e-2];%Length of Bluff Bodies [C E] in [m]

r=[1 2];%Aspect ratios
b=a./r%Minor axis half length in [m]

Mbf=rho_bf*pi*a^2*Lbf./r;%Mass of bluff body in [kg]
Jbf=Mbf.*(3.*(a./r).^2+4*Lbf^2)/12;%Mass moment of Inertia of the Bluff Bodies about the its CG [kg.m^2]
Jtip=Jbf+Mbf.*(0.5*Lbf)^2;%Mass moment of Inertia of the Bluff Bodies about the bimorph tip [kg.m^2]

rho_f = 996;%Fluid Density [kg/m^3]
M_f = Lbf*rho_f*pi.*a.^2.*b;%Fluid Added Mass [kg]

m_eq=(33/140)*Mbi+Mbf+4.*Jbf/(Lbf)^2+M_f;%Equivalent Mass [kg]
f_s=1/2/pi.*sqrt(k_eq./m_eq);%Natural Un-damped Frequency [Hz]


L_f=6e-2;%Length of Fluid in Contact
t=0:(10000-100);
V = 0.067;
v = 0.798e-3/996;
Re = 2*a*V/v; 

for i=1:2
    if (Re>=100) && (Re<1600)
        C_L=0.045+1.05*(1-Re/1600)^4.5;
    end
    if (Re>=1600) && (Re<5400)
        C_L=0.045+3*(log10(Re/1600))^4.6;
        end
     if (Re>=5400) && (Re<=225000)
        C_L=0.52-0.06*(log10(Re/1600))^-2.6;
        end
    A(i)=0.5*C_L*L_f*rho_f*2*a*V^2;
%     f_v(i)=St*V/D;
%     w_v(i)=2*pi*f_v(i);
    i=i+1;
      
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_n=sqrt(k_eq/m_eq);%Natural Un-damped Frequency [rad/s]

Lbi=76.73e-3-0.05e-3;%Bimorph Length
W= 32.98e-3-0.05e-3;%Width of the Bimorph
Hp = 0.22e-3; %Thickness of the Bimorph

Ep= 5.6e10;%Modulus of Elasticity for the PZT layers
C_p=2*epsS_33 * W * Lbi / Hp;%Effective Capacitance of the Bimorph for Parralel Connection
%k_e=sqrt((2*W*Lbi*cE_11*d_31/Lbi)^2/((2*W*Lbi*cE_11/Lbi)*(1-k_31^2)*C_p))


R_opt=2*zeta_m/k_e^2/w_n/C_p%Optimum Resistance for PEH
R_L=1e6;
%R_L = [0:1e2:(R_opt+5e4)]';
alpha=w_n.*R_L.*C_p;
w=2*pi.*[0:0.001:3]';
Omega=w./w_n;

Theta=-0.5*e_31*W*Lbi/Hp;

xx = Mbf(1)/Mbi;
cn =xx^2+0.603*xx+0.08955;
cd =xx^2+0.4637*xx+0.05718; 
mu = cn/cd;%Correction Factor for SDOF Model {3-1} Mode

num1=sqrt(1+(alpha.*Omega).^2);
deni=(1-(1+2*zeta_m.*alpha).*Omega.^2).^2;
denj=(Omega.*(2*zeta_m+(1+k_e^2).*alpha)-alpha.*Omega.^3).^2;
den1=sqrt(deni+denj);

num2=Omega;%alpha*k_e^2.*Omega;
deni=(1-(1+2*zeta_m.*alpha).*Omega.^2).^2;
denj=(Omega.*(2*zeta_m+(1+k_e^2).*alpha)-alpha.*Omega.^3).^2;
den2=sqrt(deni+denj);

num3=Omega.^2;%alpha*k_e^2.*Omega.^2;
deni=(1-(1+2*zeta_m.*alpha).*Omega.^2).^2;
denj=(Omega.*(2*zeta_m+(1+k_e^2).*alpha)-alpha.*Omega.^3).^2;
den3=deni+denj;

FFA=A(1);%Forcing Function Amplitude

%Forced Vibration Response
z=(mu/w_n^2).*(num1./den1).*FFA;%(1/k_eq).*(num1./den1).*FFA;
%Voltage Response
V=(mu*m_eq.*R_L*abs(d_31)*w_n).*(num2./den2).*FFA;%(1/abs(Theta)).*(num2./den2).*FFA;
%Power Output Response
P_out=(mu^2*m_eq.*alpha.*R_L*k_e^2./R_L./w_n).*(num3./den3).*FFA.^2;%(w_n/k_eq).*(num3./den3).*FFA.^2;
%V = abs(sqrt(P_out.*R_L));%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
plot(w/2/pi, P_out)
xlabel('Input Frequency [Hz]')
ylabel('Output Power [\muW]')
grid on
grid minor

% figure(1)
% plot(R_L, P_out)
% xlabel('Load Resistance [\Omega]')
% ylabel('Output Power [W]')
% grid on
% grid minor
% 
% figure(2)
% plot(w/2/pi, z*1000)
% xlabel('Input Frequency [Hz]')
% ylabel('Tip Displacement [mm]')
% grid on
% grid minor

figure(3)
plot(w/2/pi, V)
xlabel('Input Frequency [Hz]')
ylabel('Output Voltage across R_L = 20 k\Omega [V]')
grid on
grid minor