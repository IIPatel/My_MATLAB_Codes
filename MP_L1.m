clc
clear

Do=2e-3;
d=2e-3;%=w
fr=0.2e-3;%=to
N=20;
mu=0.75

tc=[1 2 3 4 5]';

tcavg=mean(tc)

v=pi*Do*N/60
to=fr/N
w=d
sysms a B r
r=to/tcavg
zeta=r^-1
B=atand(mu)

alph=solve(a==2*(atand(r*cosd(a)/(1-r*sind(a))))-90+B, a)
a=ans;
phi=atand(r*cosd(a)/(1-r*sind(a)))
gamma=tand(phi-a)+cotd(phi)
B
S=250e6
As=to*w/sind(phi);
Fs=S*As
Fc=(Fs*cosd(B-a))/cosd(phi+B-a)
Ff=(Fs*sind(B-a)/cosd(phi+B-a)
F=Fc*sind(a)+Ff*cosd(a)
Pc=Fc*v

