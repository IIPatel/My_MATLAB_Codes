clc
clear all
close all

Re=200;%Reynold Number
v=8.012048192771084e-07;%Kinematic Viscosity
St=[0.194 0.0636363 0.0343]';%Strohal Number for Re=200
r=[1 4 6]';%Aspect Ratios for Elliptical Cylinders (r=1;Circular)
D=2e-2.*r;%Characteristic Length of Cylinders: Minor Axis Length
f_v=St.*Re.*v./(D).^2%Vortex Shedding Frequency in Hertz for Re=200

subplot(121)
plot(r,St)
xlabel('Aspect Ratio')
ylabel('Strouhal Number')
title('r_a vs St for Re=200')
subplot(122)
plot(r, f_v)
xlabel('Aspect Ratio')
ylabel('Vortex Shedding Frequency / [Hz]')
title('r_a vs f_v for Re=200')
