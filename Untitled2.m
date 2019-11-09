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

zeta_m=0.075;%Effective Damping Ratio; Obtained from Free Vibration Tests in Fluid Medium
m_eq=0.041989170053617;%0.580630807575161;%Equivalent Mass of the Harvester
k_eq=3.048593557972514e+02;%Equivalent Stiffness of the Harvester
w_n=sqrt(k_eq/m_eq);%Natural Un-damped Frequency [rad/s]

Lbi=76.73e-3-0.05e-3;%Bimorph Length
W= 32.98e-3-0.05e-3;%Width of the Bimorph
Hp = 0.22e-3; %Thickness of the Bimorph

Ep= 5.6e10;%Modulus of Elasticity for the PZT layers
C_p=2*epsS_33 * W * Lbi / Hp;%Effective Capacitance of the Bimorph for Parralel Connection
%k_e=sqrt((2*W*Lbi*cE_11*d_31/Lbi)^2/((2*W*Lbi*cE_11/Lbi)*(1-k_31^2)*C_p))


R_opt=2*zeta_m/k_e^2/w_n/C_p;%Optimum Resistance for PEH
R_L=R_opt

alpha=w_n*R_L*C_p;
w=2*pi.*[0:0.001:20]';
Omega=w./w_n;

I = 7.8894e-13;
Theta= -d_31*k_eq*W*Lbi*(Lbi+10e-2)*Hp/I;

num1=sqrt(1+(alpha.*Omega).^2);
deni=(1-(1+2*zeta_m*alpha).*Omega.^2).^2;
denj=(Omega.*(2*zeta_m+(1+k_e^2)*alpha)-alpha.*Omega.^3).^2;
den1=sqrt(deni+denj);

num2=alpha*k_e^2.*Omega;
deni=(1-(1+2*zeta_m*alpha).*Omega.^2).^2;
denj=(Omega.*(2*zeta_m+(1+k_e^2)*alpha)-alpha.*Omega.^3).^2;
den2=sqrt(deni+denj);

num3=alpha*k_e^2.*Omega.^2;
deni=(1-(1+2*zeta_m*alpha).*Omega.^2).^2;
denj=(Omega.*(2*zeta_m+(1+k_e^2)*alpha)-alpha.*Omega.^3).^2;
den3=deni+denj;

FFA=0.6;%Forcing Function Amplitude

%Forced Vibration Response
z=(1/k_eq).*(num1./den1).*FFA;
%Voltage Response
V=(1/abs(Theta)).*(num2./den2).*FFA;
%Power Output Response
P_out=(w_n/k_eq).*(num3./den3).*FFA^2;
V = sqrt(P_out.*R_L);

figure(1)
plot(w/2/pi, P_out*1000)
xlabel('Input Frequency [Hz]')
ylabel('Output Power [mW]')
grid on
grid minor

figure(2)
plot(w/2/pi, z*1000)
xlabel('Input Frequency [Hz]')
ylabel('Tip Displacement [mm]')
grid on
grid minor

figure(3)
plot(w/2/pi, V)
xlabel('Input Frequency [Hz]')
ylabel('Output Voltage across R_L =  k\Omega [V]')
grid on
grid minor