clc
clear all
close all

num0=0.1;
den0=[1 1 0];
Gp=tf(num0,den0);% Vehicle T.F.

num1=10;
den1=[1 5];
Ga=tf(num1,den1);% Actuator T.F.

num2=100;
den2=[1 20 100];
H=tf(num2,den2);% Accelerometer T.F.

Kp=6.23166444170877;% Proportional gain
Ki=0;% Integral gain
Kd=10.9503477881035;% Derivative gain

num3=[Kd Kp Ki];
den3=[1 0];
PID=tf(num3,den3)% PID T.F.

S1D=series(H, PID);
S2D=series(S1D, Ga);
FD=feedback(Gp,S2D);
XD=-FD;%X(s)/D(s)

S1R=series(Ga, PID);
S2D=series(S1R, Gp);
FR=feedback(S2D,H);
XR=FR;%X(s)/R(s)

R=tf(1,[1 0]);%Unit Step R(s)
D=R;%D(s)
SP=XR*R+XD*D;%Superposition to obtainX(s)
E=R-SP;%E(s)
U=E*Ga*PID;%U(s)

step(XD); hold on; step(XR); ; impulse(SP); impulse(U)
legend('R(s)=0','D(s)=0','Combined','U(s)')