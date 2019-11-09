clc; clear all; close all;

y_0(1)=0.5;
y(1)=0.5;
h=0.5;
t=[0:h:2]';

for i=1:(length(t)-1)
fxny_i=diffeq_E2(t(i),y(i))
y_0(i+1)=y(i)+fxny_i*h
y_dash=diffeq_E2(t(i+1),y_0(i+1))
y_bari=(fxny_i+y_dash)/2
y(i+1)=y(i)+y_bari*h
end

