clc
clear all
close all

n = 'A0057CH1.csv'% Name of Sample File
A = xlsread(n,'A17:A4017');
xcf = xlsread(n,'B9:B9');%ms/div
ycf = xlsread(n,'B6:B6');%V/div
pcf = 25;%points/div

A = A.*(ycf/pcf);%Scaling points according to the scale
B = A;
t = 0:length(A)-1;
t = t'.*(xcf/pcf);
     
V1 = [10.6 15 15 14];
V2 = [4.8 5.8 6.6 5.8];
i = 1;

 while i <= length(V1)
zeta(i) = log(V1(i)/V2(i))/sqrt(4*pi+(log(V1(i)/V2(i)))^2)
     i = i+1;
 end

 zeta = sum(zeta)/length(zeta)
 
figure(1)
plot(t,A(:,1))
title('Experimental Data')
xlabel('Time [s]')
ylabel('Voltage [V]')
grid on
grid minor



