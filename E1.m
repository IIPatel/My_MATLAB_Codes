clc; clear all; close all;

y_0(1)=2;
y(1)=2;
h=1;
x=[0:4]';

for i=1:(length(x)-1)
fxny_i=diffeq_E1(x(i),y(i))
y_0(i+1)=y(i)+fxny_i*h
y_dash=diffeq_E1(x(i+1),y_0(i+1))
y_bari=(fxny_i+y_dash)/2
y(i+1)=y(i)+y_bari*h
end

