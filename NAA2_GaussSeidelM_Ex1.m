clc; clear all; close all;
 
b=[27 -61.5 -21.5]';
x=[0 0 0]';

 
xi=[0 0 0]';
er_x=[0 0 0]';

it=0;
 
while it<5
x(1)=(b(1)-2*x(2)+x(3))/10;
 
x(2)=(b(2)+3*x(1)-2*x(3))/(-6);
 
x(3)=(b(3)-x(2)-x(1))/5;
 
er_x(1)=abs((xi(1)-x(1))/x(1));
er_x(2)=abs((xi(2)-x(2))/x(2));
er_x(3)=abs((xi(3)-x(3))/x(3));

xi=x;
x;
er_x=er_x*100;
table=[x(1) er_x(1) x(2) er_x(2) x(3) er_x(3)]
it=it+1;
end
