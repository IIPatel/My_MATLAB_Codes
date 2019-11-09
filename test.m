clc
clear all
close all
d_31= -186e-12;
d_33=6.70e-010;
sE_33=2.07e-011;
sE_11=1.65e-011;
sE_12=-4.78e-012;
sE_13=-8.45e-012;
sE_31=-8.45e-012;
cE_11=sE_11/(sE_11^2-sE_12^2);
cE_33=sE_33^-1 
k_p=0.65;%material data
sigmaP = -sE_12/sE_11;%Planar Poisson Ratio
epsT_33=8.854187817e-12*3400;

epsS_33 = epsT_33 - 2*d_31^2/(sE_11+sE_12);%Permitivity at constant strain
%epsS_33=epsS_33/8.854187817e-12;
e_31=d_31/(sE_11+sE_12)

%e_33=d_31/sE_11;
%epsS_33 = epsT_33 - d_31*e_33;%Permitivity at constant strain

k_31=sqrt(k_p^2*0.5*(1-sigmaP))
k_e=sqrt(k_31^2./(1-k_31^2))
k_e=sqrt((e_31^2./(cE_33*epsS_33))/(1-(e_31^2./(cE_33*epsS_33))))

%e_31=sqrt(k_e^2*sE_33^-1*epsS_33)

zeta_m=0.075;
m_eq=0.580630807575161;
k_eq=3.048593557972514e+02;
w_n=sqrt(k_eq/m_eq);

Lbi=76.73e-3-0.05e-3;
W= 32.98e-3-0.05e-3;%Width of the piezobeam
Hp = 0.22e-3; %Thickness of the piezo

Ep= 5.6e10;
C_p=2*epsS_33 * W * Lbi / Hp
%k_e=sqrt((2*W*Lbi*cE_11*d_31/Lbi)^2/((2*W*Lbi*cE_11/Lbi)*(1-k_31^2)*C_p))


R_opt=2*zeta_m/w_n/C_p/(sqrt(4*zeta_m^2+k_31^4))
R_L=R_opt

alpha=w_n*R_L*C_p;
w=2*pi.*[0:0.001:7]';
Omega=w./w_n;

Theta=-e_31*2*W*Lbi/Hp;

num1=sqrt(1+(alpha.*Omega).^2);
deni=(1-(1+2*zeta_m*alpha).*Omega.^2).^2;
denj=(Omega.*(2*zeta_m+(1+k_e^2)*alpha)-alpha.*Omega.^3).^2;
den1=sqrt(deni+denj);

k2=(Lbi^2/(3*W))*(2*Lbi+1.5*10e-2)/(2*Lbi+10e-2-Lbi);
num2=j.*w.*(Ep*d_31*Hp/epsT_33/k2);
deni=w_n^2/R_L/C_p-((R_L*C_p)^-1+2*zeta_m.*w_n).*w.^2;
denj=j.*w.*(w_n^2*(1+k_31^2)+2*zeta_m*w_n/R_L/C_p-w.^2);
den2=(deni+denj);

num3=alpha*k_e^2.*Omega.^2;
deni=(1-(1+2*zeta_m*alpha).*Omega.^2).^2;
denj=(Omega.*(2*zeta_m+(1+k_e^2)*alpha)-alpha.*Omega.^3).^2;
den3=deni+denj;

FFA=0.6;%Forcing Function Amplitude

%Forced Vibration Response
z=(1/k_eq).*(num1./den1).*FFA;
%Voltage Response
V=(num2./den2).*FFA;
%Power Output Response
P_out=(abs(V)).^2./R_L;

subplot(221)
plot(w/2/pi, P_out)
xlabel('Input Frequency [Hz]')
ylabel('Output Power [W]')
subplot(222)
plot(w/2/pi, z*100)
xlabel('Input Frequency [Hz]')
ylabel('Tip Displacement [cm]')
subplot(223)
plot(w/2/pi, abs(V))
xlabel('Input Frequency [Hz]')
ylabel('Output Voltage across R_L = 55 k\Omega [V]')
