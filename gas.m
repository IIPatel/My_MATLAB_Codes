clc
clear

P2_P1=2.25;

%Rho2_Rho1=4.5/1.5;

syms M1 k

k=1.4

Mi= solve(P2_P1==(2*k*M1^2-(k-1))/(k+1), M1)

Mi=406^(1/2)/14

a1=sqrt(k*287*308)

Mi*a1

T2_T1=((2*k*Mi^2-(k-1))*(2+(k-1)*Mi^2))/(((k+1)^2)*(Mi^2))

Mf=sqrt(((k-1)*Mi^2+2)/(2*k*Mi^2-(k-1)))

T2=T2_T1*(35+273);

a2=sqrt(k*287*T2);

V2=Mi*sqrt(k*287*308)-Mf*a2;

Mi=1.3817;

Mf=sqrt(((k-1)*Mi^2+2)/(2*k*Mi^2-(k-1)))

T2_T1=((2*k*Mi^2-(k-1))*(2+(k-1)*Mi^2))/(((k+1)^2)*(Mi^2));

M=(V2/a2)+Mf*sqrt(T2_T1)

((M*a2-V2)/Mf)^2/(k*287)

T2_T1*(T2)

P2_P1=(2*k*M^2-(k-1))/(k+1)

P2_P1*270
