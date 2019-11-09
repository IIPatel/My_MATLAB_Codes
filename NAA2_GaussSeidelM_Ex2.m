clc; clear all; close all;
 
b=[3800 1200 2350]';
c=[0 0 0]';

 
ci=[0 0 0]';
er_c=[0 0 0]';

it=0;
 
while it<4
c(1)=(b(1)+3*c(2)+c(3))/15;
 
c(2)=(b(2)+3*c(1)+6*c(3))/18;
 
c(3)=(b(3)+4*c(1)+c(2))/12;
 
er_c(1)=abs((ci(1)-c(1))/c(1));
er_c(2)=abs((ci(2)-c(2))/c(2));
er_c(3)=abs((ci(3)-c(3))/c(3));

ci=c;
c;
er_c=er_c*100;
table=[c(1) er_c(1) c(2) er_c(2) c(3) er_c(3)]
it=it+1;
end
