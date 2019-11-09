clc
clear all 
close all

t = [0:2*pi/1000:2*pi];
x =5+cos(5.*t).*cos(t);
y=5+cos(5.*t).*sin(t);

subplot(231)
plot(x,y)
xlabel('X-coordinate')
ylabel('Y-coordinate')
title('5-petal flower')
subplot(232)
plot(t,x)
xlabel('Phase')
ylabel('X-coordinate')
subplot(233)
plot(t,y)
xlabel('Phase')
ylabel('Y-ccordinate')

%Fourier Series Approximation
%Omega
wx=2;
%Coefficients (obtained from the Curve Fitting Tool)
wx=2;
ax=[10 0 0.5 0.5];
bx=[0 0 0];

X=ax(1)/2+ax(2)*cos(1*wx.*t)+ax(3)*cos(2*wx.*t)+ax(4)*cos(3*wx.*t);
subplot(234)
plot(t,X)
xlabel('Phase')
ylabel('X-coordinate')
title('n=3')

wy=2;
ay=[10 0 0 0];
by=[0 -0.5 0.5];

Y=ay(1)/2+by(1)*sin(1*wx.*t)+by(2)*sin(2*wx.*t)+by(3)*sin(3*wx.*t);

subplot(235)
plot(t,Y)
xlabel('Phase')
ylabel('Y-coordinate')
title('n=3')
subplot(236)
plot(X,Y)
xlabel('X-coordinate')
ylabel('Y-coordinate')
title('5-petal flower, n=3')

X2=ax(1)/2+ax(2)*cos(1*wx.*t)+ax(3)*cos(2*wx.*t);
Y2=ay(1)/2+by(1)*sin(1*wx.*t)+by(2)*sin(2*wx.*t);

X1=ax(1)/2+ax(2)*cos(1*wx.*t);
Y1=ay(1)/2+by(1)*sin(1*wx.*t);
subplot(231)
hold on
plot(X2,Y2)
hold on
plot(X1,Y1)
legend('Original', 'n=2', 'n=1')

