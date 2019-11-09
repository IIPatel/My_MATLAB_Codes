clc; clear all; close all;
 
t=[5 8 12]';
v=[106.8 177.2 279.2]';
a=[1 2 5]';
 
ai=[1 2 5]';
era=[0 0 0]';
 
it=0;
 
while it<2
a(1)=(v(1)-a(3)-(a(2)*t(1)))/(t(1)^2);
 
a(2)=(v(2)-a(3)-(a(1)*t(2)^2))/t(2);
 
a(3)=v(3)-a(1)*t(3)^2-a(2)*t(3);
 
era(1)=abs((ai(1)-a(1))/a(1));
era(2)=abs((ai(2)-a(2))/a(2));
era(3)=abs((ai(3)-a(3))/a(3));
ai=a;
a
era=era*100
it=it+1;
end
