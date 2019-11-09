clc
clear all
close all
%PVEH system using SDOF Stress-Voltage Analogous model
Lbi=76.73e-3-0.05e-3; %Length of the overhung piezobeam
Lp=Lbi-1.5e-2;
W= 32.98e-3-0.05e-3;%Width of the piezobeam
Hs= 0.22e-3;%Thickness of the shim
Hp = 0.22e-3; %Thickness of the piezo
Es = 11.2e10; sm = 1/Es; %Young’s Modulus and compliant of the shim
Ep = 5.6e10; s11 = 1/Ep; %Young’s Modulus and compliant of piezo
rho_s = 8780; %Density of shim
rho_p = 7500; %Density of piezo
d31 = -186e-12; %-190e-12, Piezoconstant;
zeta = 0.008; %damping ratio
yc = (Hp+Hs)/2; %Location of neutral axis
mass_beampiezo = 2*rho_p*Lp*W*Hp; 
mass_beamshim = rho_s*Lbi*W*Hs;
Mbi = mass_beampiezo + mass_beamshim%Total mass of EH beam
mass_dist = Mbi/Lbi; % Distributed mass per length
permit_cons_elect = 8.854187817e-12 * 3400; %Permitivity in free space
perm_cons_strain = (permit_cons_elect - d31^2*Ep);
%perm_cons_strain = perm_cons_strain/8.854187817e-12%Permitivity at constantstrain
R = 1e-3; % Resistance value in Ohms

Ip=2*W*((Hp^3)/12+Hp*(0.5*Hp+0.5*Hs)^2);%Second moment of area about the neutral axis of the piezo layers.
Is=(W*Hs^3)/12;%Second moment of area about the neutral axis of the shim.
EI=Ep*Ip+Es*Is;%Effective Flexural Stiffness of the bimorph

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_eq=3*EI/Lbi^3;%Equivalent Spring Stiffness of the solid structure

%Bluff body parameters
rho_bf=7850;%Bluff Body material density / kg m^-3
a=1e-2;%Major axis half length in m
Lbf=10e-2;%Length of Bluff Body in m
Le=Lbf;%Electric Field Length (Piezo Layer Length)
r=[1];%Aspect ratios
b=a./r.*100%Minor axis half length in cm
Mbf=rho_bf*pi*a^2*Lbf./r%Mass of bluff body in kg
Jbf=Mbf.*(3.*(a./r).^2+4*Lbf^2)/12;%Mass moment of Inertia of the Bluff Body about the its CG
Jtip=Jbf+Mbf.*(0.5*Lbf)^2;%Mass moment of Inertia of the Bluff Body about the bimorph tip ;

m_eq=(33/140)*Mbi%+Mbf+4.*Jbf/(Lbf)^2;%Equivalent Mass

rho=1000;%Density of water kg m^-3
v=0.798e-3/rho;%Kinematic Viscosity of Water @ 30 degC
D=2e-2;%Circular Cylinder Diameter/m
St=0.22;%St for circular cylinder 

L_f=15e7;%Length of Fluid in Contact
t=0:(10000-100);


for i=1:(10000-99)
    
    Re=99+i;
    V=Re*v/D;
    
    if (Re>=100) && (Re<1600)
        C_L=0.045+1.05*(1-Re/1600)^4.5;
    end
    if (Re>=1600) && (Re<5400)
        C_L=0.045+3*(log10(Re/1600))^4.6;
        end
     if (Re>=5400) && (Re<=225000)
        C_L=0.52-0.06*(log10(Re/1600))^-2.6;
        end
    A(i)=0.5*C_L*L_f*rho*D*V^2;
    f_v(i)=St*V/D;
    w_v(i)=2*pi*f_v(i);
    i=i+1;
      
end
        
T_v=1./f_v;
Lift=A.*cos(w_v.*t);
w=w_v; %Excitation frequency range in Radians
w_hz = w./2/pi; %Excitation frequency in Hz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

% Mechanical modal constants for first mode of clamp-free beam, 
%Equations (3.8), (3.9)
Lamda = 1.87510407; % Standard value for cantilever beam (Mode shape books)
Sigma = 0.734095514; % % Standard value for cantilever beam
phi_r = ((cosh(Lamda)-cos(Lamda))-Sigma*(sinh(Lamda)-sin(Lamda)))/sqrt(Mbi); 
%Part of Equation %(3.8), terms within brackets
trans_constant_1 = 2*Sigma/Lamda; %%%Gamma translation
trans_constant_2 = sqrt(Lbi/mass_dist);
trans_constant = trans_constant_1*trans_constant_2;
phi_deriv1 = (sinh(Lamda)+ sin(Lamda))-Sigma*(cosh(Lamda)-cos(Lamda));
phi_derivative = Lamda*phi_deriv1/(sqrt(Mbi)*Lbi)

% Electromechanical constants in mechanical domain %

elect_const = -(Ep*d31*W*(Hp + Hs))/(2); % Equation (3.16)
Xr = elect_const * phi_derivative; %Equation (3.18)

I=(2*W*Hp^3/12+W*2*Hp*Hp^2)+(Ep/Es)*W*Hs^3/12;
b_aa=2*I/(Hp*(2*Lbi+Lbf-Le));
b_a=(3*Hp/Lbi^2)*((2*Lbi+Lbf-Le)/(2*Lbi+1.5*Lbf));

% Calculating natural frequencies of the PV energy harvester
w_r = sqrt(k_eq/m_eq);%Lamda^2*sqrt(EI/(mass_beam*Lbi^3)); %Equivalence of equation (3.11)
w_r_Hz = w_r/(2*pi) %Natural frequency in Hz

% Voltage constant, V and Xr, electromechanical constants in mechanical domain
% Capacitance of piezoelectric beam
Cp =2*perm_cons_strain * W * Lbi / (Hp); %permit_cons_elect*2 * W * Lbi / (Hp); % Equation (3.24), 2 for series connection
modal_const = -d31*Ep*yc*W*phi_derivative; %
K=sqrt(Ep*d31/permit_cons_elect);
K=abs(K);

%% Term-by-term voltage FRF calculation
Volt_nume =-j*w.*Ep*d31*Hp*b_a/permit_cons_elect;
Volt_denom1= ((R*Cp)^-1)*w_r^2-((R*Cp)^-1+2*zeta*w_r).*w.^2;
Volt_denom2 = j*w.*(w_r^2*(1+K^2)+2*zeta*w_r/(R*Cp)-w.^2);
Volt_denom_final = (Volt_denom2 + Volt_denom1);

% Voltage FRF
VOLTAGE_FRF = (Volt_nume./Volt_denom_final.*A);
VOLTAGE_FRF_abs = abs(VOLTAGE_FRF);
plot(w_hz,VOLTAGE_FRF_abs,'k')
%semilogy(w_hz,VOLTAGE_FRF,'r')% for semilog “y” axis
title('R = 150 [k\Omega], Voltage')
xlabel('Frequency / (Hz)')
ylabel('Volt')
%axis([0 200 0 2])% Specifying the range on the axis if needed
figure;
hold on

% Calculating CURRENT (mA)FRFs for PVEH System %
Current_FRF = VOLTAGE_FRF_abs./R*1000;
title('CURRENT FRFs')
xlabel('Frequency Hz')
ylabel('Current / (mA)')
semilogy(w_hz,Current_FRF);
figure;
hold on

%Calculating POWER FRFs for PVEH system normalised by g of acceleration
Power = (VOLTAGE_FRF_abs.^2./R)*1000; % milliWatts default was comlex Voltage
plot(w_hz, Power)
title('Power')
xlabel('Frequency / (Hz)')
ylabel('mW')
% It is important to note that the program is valid for a particular case based on equations of Chap. 3. It can % be different for different piezo materials, input
%frequencies, resonance frequencies, connected % resistor % value, damping
%value and other piezoelectric constants.