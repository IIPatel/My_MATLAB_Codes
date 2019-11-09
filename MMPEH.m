clc
clear all
close all
%PVEH system using distributed parameter model
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
w=0:0.001:2000; %Excitation frequency range in Radians
w_hz = w./2/pi; %Excitation frequency in Hz
mass_beampiezo = 2*rho_p*Lp*W*Hp; 
mass_beamshim = rho_s*Lbi*W*Hs;
mass_beam = mass_beampiezo + mass_beamshim; %Total mass of EH beam
mass_dist = mass_beam/Lbi; % Distributed mass per length
permit_cons_elect = 8.854187817e-12 * 3400; %Permitivity in free space
perm_cons_strain = (permit_cons_elect - d31^2*Ep); %Permitivity at constantstrain
R = 12e3; % Resistance value in Ohms

% Mechanical modal constants for first mode of clamp-free beam, 
%Equations (3.8), (3.9)
Lamda = 1.87510407; % Standard value for cantilever beam (Mode shape books)
Sigma = 0.734095514; % % Standard value for cantilever beam
phi_r = ((cosh(Lamda)-cos(Lamda))-Sigma*(sinh(Lamda)-sin(Lamda)))/sqrt(mass_beam); 
%Part of Equation %(3.8), terms within brackets
trans_constant_1 = 2*Sigma/Lamda; %%%Gamma translation
trans_constant_2 = sqrt(Lbi/mass_dist);
trans_constant = trans_constant_1*trans_constant_2;
phi_deriv1 = (sinh(Lamda)+ sin(Lamda))-Sigma*(cosh(Lamda)-cos(Lamda));
phi_derivative = Lamda*phi_deriv1/(sqrt(mass_beam)*Lbi);

% Electromechanical constants in mechanical domain %

Ip=2*W*((Hp^3)/12+Hp*(0.5*Hp+0.5*Hs)^2);%Second moment of area about the neutral axis of the piezo layers.
Is=(W*Hs^3)/12;%Second moment of area about the neutral axis of the shim.
EI=Ep*Ip+Es*Is;%Effective Flexural Stiffness of the bimorph
elect_const = -(Ep*d31*W*(Hp + Hs))/(2); % Equation (3.16)
Xr = elect_const * phi_derivative; %Equation (3.18)

% Calculating natural frequencies of the PV energy harvester
w_r = Lamda^2*sqrt(EI/(mass_beam*Lbi^3)); %Equivalence of equation (3.11)
w_r_Hz = w_r/2/pi %Natural frequency in Hz

% Voltage constant, V and Xr, electromechanical constants in mechanical domain
% Capacitance of piezoelectric beam
Cp = perm_cons_strain * W * Lbi / (Hp); % Equation (3.24), 2 for series connection
modal_const = -d31*Ep*yc*W*phi_derivative; %
%% Term-by-term voltage FRF calculation for single-mode expression, Equation (3.31)
Volt_denom1 = (w_r^2 - w.^2) + (j*2*w.*w_r*zeta);
Volt_nume =-j*2*w.*R*mass_dist*trans_constant*modal_const;
Volt_denom11 = (j*2*w.*R*Xr*modal_const);
Volt_denom22 = (2+j*w.*R*Cp);
Volt_denom33 = Volt_denom1.*Volt_denom22;
Volt_denom44 = Volt_denom22.*Volt_denom1;
Volt_denom_final = (Volt_denom44 + Volt_denom11);

% Voltage FRF as per equation (3.31), for per “g”, multiply by 9.81
VOLTAGE_FRF = (Volt_nume./Volt_denom_final);
VOLTAGE_FRF_abs = abs(VOLTAGE_FRF);
plot(w_hz,VOLTAGE_FRF_abs,'k')
%semilogy(w_hz,VOLTAGE_FRF,'r')% for semilog “y” axis
title('R = 12 [k\Omega], Voltage')
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