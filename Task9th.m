clc; clear all; close all;
x=[0 0 0 0]';
t=[2 1 0 -2]';
error=1;
it=0;

while error > 10e-3
 
x1=(6+2*x(2)+x(3)-x(4))/5;
x2=(2*x(1)-x(3))/4;
x3=(6-x(1)-2*x(2)+x(4))/6;
x4=(-14+x(1)-x(3))/6;

xi_plus1=[x1 x2 x3 x4]';
error=abs((xi_plus1-x)./x);
x=[x1 x2 x3 x4]';
it=it+1;
end
fprintf('The number of iterations required are: %d', it)

