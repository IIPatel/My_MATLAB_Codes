clc
clear all
close all

Li=127e-3;
L0=130e-3;
k=0.17e3;
Do=12.7e-3;Di=11.48e-3;
d=1.22e-3;
N=85;
g=9.81;

mc=600e-3;
ms=28.6e-3;
md=[0 0.2 0.4 0.6]; 
m=mc+md;
m=m';
meff=ms./3;

%S1=g./k;
%L=S1.*m+L0;

M=mc+meff+md;
wn=sqrt(k./M);
t=2*pi./wn;

%S2=4*pi^2./k;
%Gamma=S2.*(mc+meff);
%tsq=S2.*m+Gamma;

G=80e9;
D=(Di+Do)/2;
k3=d^4*G/(8*D^3*N)
y=[0 5 1015 20]';
L=y+L0;
ki1=(m.*g)./y

subplot(211)
f1=fit(m,L,'poly1');
plot(f1,m,L);
xlabel('Total Mass of Disk(s) / kg')
ylabel('Extension of Spring / mm')
title('Load-Extension Graph')

tsq=[0 5 10 20]';
ki2=4*pi^2./tsq
subplot(212)
f2=fit(m,tsq,'poly1');
plot(f2,m,tsq);
xlabel('Total Mass of Disk(s) / kg')
ylabel('\tau^2 / s^2')
title('Oscillation Period^2 against Load Graph')

%k1=g.*m./S1
%k2=4*pi^2./S2

% k3-k1
% k3-k2
% k2-k1

