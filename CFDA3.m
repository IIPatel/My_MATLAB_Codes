clc
clear all
close all

A=[2125 -125 0 0 0
   -125 2250 -125 0 0
   0 -125 2250 -125 0
   0 0 -125 2250 -125
   0 0 0 -125 2375];%Coefficient Matrix for Ax=B

To=[199.999569147615	199.992675509460	199.868590022671	197.641944898624	157.686418152559]';%Update with the row vector of x from Workspace
B=[1875*To(1)+125*To(2)
    125*To(1)+1750*To(2)+125*To(3)
    125*To(2)+1750*To(3)+125*To(4)
    125*To(3)+1750*To(4)+125*To(5)
    124*To(4)+1625*To(5)]%Column Vector

x=A\B;%Solve Linear System
x=x';%Convert to row vector form

L=[0 0.002 0.006 0.010 0.014 0.018 0.02]';
T=[200 x 0]';
    
plot(L, T)
title('Temperature Distribution across Thin Plate (L = 2 cm) at t = 2 s')
xlabel('Length [m]')
ylabel('Temperature [^oC]')
grid on
grid minor