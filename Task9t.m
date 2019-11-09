clc; clear all; close all;
 
v=[106.8 177.2 279.2]';
a=[0 0 0]';

 
ai=[0 0 0]';
er_a=[0 0 0]';
er_s=[5 5 5]';

it=0;
 
while er_a<=er_s
a(1)=(v(1)-a(3)-(a(2)*t(1)))/(t(1)^2);
 
a(2)=(v(2)-a(3)-(a(1)*t(2)^2))/t(2);
 
a(3)=v(3)-a(1)*t(3)^2-a(2)*t(3);
 
er_a(1)=abs((ai(1)-a(1))/a(1));
er_a(2)=abs((ai(2)-a(2))/a(2));
er_a(3)=abs((ai(3)-a(3))/a(3));

ai=a;
a
er_a=er_a*100
it=it+1;
end
