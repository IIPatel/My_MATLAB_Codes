clc; clear all; close all;
R=[0:15000];
T=(1.129241e-3+2.341077e-4.*log(R)+8.775468e-8.*(log(R)).^3).^-1;
plot(R, T )
grid on;
xlabel('Resistance/Ohms')
ylabel('Temperature/Kelvin')
title('10K3A Betatherm Thermistor T-R Charecteristic Curve')

xi=14000;
it=1;

while it<=3
    [resi_xi, dresi_xi]=resi(xi);
    xi_plus1=xi-(resi_xi./dresi_xi);
    err=abs((xi_plus1-xi)/(xi_plus1))*100;
    xi=xi_plus1
    err
    it=it+1;
end
