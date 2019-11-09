clc
clear

format long

E=210e9;
b=19.75e-3;
h=4.74e-3;

L=[5 6 7 9 15 10 11 12 13]';
mV=[0.138 0.171 0.205 0.271 0.476 0.322 0.356 0.388 0.419]';
mV=mV./1000;

eps=mV.*2.05^-1
eps_th=12.*(L.*0.25*h)./(2*b*h^3*E)

sigma=eps.*E
sigma_th=L.*0.25/(74.6e-9)

plot(eps, sigma)
%plot(eps, sigma, eps_th, sigma_th)