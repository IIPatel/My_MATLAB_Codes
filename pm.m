clc
clear all 
close all

%%I=[0.0018:0.000001:0.00185];
%P=300.*((1-((1+I).^(-120)))./I)
%%plot(I,P)
%%xlabel('i/mo^-1')
%%ylabel('PV/$')

%Q_v=[2180 2251.6 2340]';
%Q_k=[168 169.4 171]';

plot(Q_k, Q_v)
xlabel('Cooling Capacity due to Convection/W')
ylabel('Cooling Capacity due to Evaporation/W')
title('Cooling Capacity due to Evaporation vs Convection')

