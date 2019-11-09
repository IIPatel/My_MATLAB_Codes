clc
clear all 
close all
format long

ms = 120;
P = [12.5 10 7 5 2 1 0.5 0.1 0.05];%Extraction Pressure
h2 = [186.405 183.895 180.879 178.866 175.844 174.835 174.331 173.927 173.877];
h3 = [1511.46 1407.87 1267.44 1154.5 908.622 762.683 640.185 417.436 340.476];
h6 = [3383.69 3311.93 3204.29 3110.22 2887.09 2746.04 2618.82 2360.92 2263.6];
cp = 4.18;
dT = 18;
i = 1;
while i <= length(P)
%Specific Enthalpies at respective states [kJ/kg]

h = [173.852
h2(i)%h(2) to vary
h3(i)%h(3) to vary
779.508
3467.28
h6(i)%h(6) to vary
2037.84]';

y(i) = (h(3)-h(2))/(h(6)-h(2));
Wnet(i) = ms*(y(i)*(h(5)-h(6))+(1-y(i))*((h(5)-h(7))-(h(2)-h(1)))-(h(4)-h(3)));
Qb(i) = ms*(h(5)-h(4))
eta(i) = (Wnet(i)/Qb(i))*100
mcw(i) = ms*((1-y(i))*(h(7)-h(1)))/(cp*dT)

i = i+1;
end 

figure(1)
plot(P,eta)
title('Cycle Performance Variation with OFWH Operating Pressure')
xlabel('OFWH Pressure [MPa]')
ylabel('Cycle Efficieny (\eta_t_h) [%]')
grid on
grid minor

%Pressure Variation
eta = 0.8;
P = [8.5 7 5 3 2 1 0.7 0.5 0.1 0.05];%Reheat Pressure
h4s = [3272.4 3214.39 3119.54 2988.49 2894.65 2752.65 2685.63 2625.02 2366.36 2268.77];
h5 = [3342.96 3362.48 3387.71 3412.1 3424.01 3435.74 3439.22 3441.54 3446.15 3446.73];
h6s = [2039.99 2073.83 2129.73 2210.21 2271.76 2374.31 2426.3 2475.04 2724.16 2860.99];
i = 1;
while i <= length(P)
% Specific Enthalpies at respective states [kJ/kg]

h = [151.494
161.518%2s
0
3322.89
h4s(i)%h(5) to vary
0%h(6) to vary
h5(i)%h(7) to vary
h6s(i)%6s
0]';

h(3) = (h(2)-h(1))/eta+h(1);
h(6) = h(4)-eta*(h(4)-h(5));
h(9) = h(7)-eta*(h(7)-h(8));

qin(i) = (h(4)-h(3))+(h(7)-h(6));
qw(i) = (h(9)-h(1));
etaP(i) = (1-qw(i)/qin(i))*100


i = i+1;
end 

figure(3)
plot(P,etaP)
title('Cycle Performance Variation with Reheat Pressure')
xlabel('Reheat Pressure [MPa]')
ylabel('Cycle Efficieny (\eta_t_h) [%]')
grid on
grid minor

%Temperature Variation for P = 0.7 MPa 
eta = 0.8;
T = [170 180 200 250 300 350 400 480 550 600];%Reheat Temperature
h5 = [2775.36 2799.38 2845.29 2954.12 3059.5 3164.13 3269.14 3439.22 3590.82 3700.9];
h6s = [2073.77 2090.35 2121.02 2188.68 2248.2 2302.33 2352.47 2426.3 2485.82 2525.98];
i = 1;
while i <= length(T)
% Specific Enthalpies at respective states [kJ/kg]
h = [151.494
161.518%2s
0
3322.89
2685.63
0
h5(i)%h(7) to vary
h6s(i)%6s varies as a result
0]';

h(3) = (h(2)-h(1))/eta+h(1);
h(6) = h(4)-eta*(h(4)-h(5));
h(9) = h(7)-eta*(h(7)-h(8));

qin(i) = (h(4)-h(3))+(h(7)-h(6));
qw(i) = (h(9)-h(1));
etaT(i) = (1-qw(i)/qin(i))*100


i = i+1;
end 

figure(4)
plot(T,etaT)
title('Cycle Performance Variation with Reheat Exit Temperature for P = 0.7 MPa')
xlabel('Reheat Exit Temperature (T_5) [^oC]')
ylabel('Cycle Efficieny (\eta_t_h) [%]')
grid on
grid minor
