clc
clear
format long

a=0.116380330080806;
a=a*pi/180
Do=34e-3;
d=0.5e-3;
f=0.111e-3;
N=850;
mu=[0.4:-0.1:0.1]'

zeta=exp(mu.*pi/2-mu.*a)

r=1./zeta;

phi=atan(r.*cos(a)./(1-r.*sin(a)))
phi=phi.*180/pi

gamma=tand(phi-a)+cotd(phi)

v=pi*Do*N/60

RMR=v*f*d

tc=zeta.*f