clc; clear all; close all;

%Given Application Parameters:
k=106663.56011183596; %Spring Stiffness (in N/m)
Fmax=5886;% MaximumForce (in N)
OD=0.118;% Coil Outer Diameter (in m)
L0=0.2592;%Free Length(in m)
Xi=0.05;%Spring Overrun to Closure Ratio [Robust Linearity]
Fmin=981;%Minimum Force (in N)
Fa=(Fmax-Fmin)/2;%Amplitude of Force (in N)
Fm=(Fmax+Fmin)/2;%Mean Force (in N)
Ssa=398e6;%from Zimmereli's Data [Peened] (in Pa)
Ssm=534e6;%from Zimmereli's Data [Peened] (in Pa)

%Material Properties and/or Parameters:
E=196.5e9;%Modulus of Elasticity (in Pa)
G=78.6e9;%Shear Modulus (in Pa)
gamma=76.518e3;%Specific Weight of Wire Material (in N/m^3)
RC=1.0;%Relative Cost of Wire Material

%Data Conversion & Computation:
d=[0.250
0.312
0.375
0.438
0.500
0.562
0.625];%Wire diameters (in in)

d=d.*0.0254;%Wire diameter (in m)

D=OD-d;% Mean Coil diameter (in m)

C=D./d;%Spring Index

KB=(4.*C+2)./(4.*C-3);% Stress Correction Factor

Taumax=KB.*((8.*Fmax.*D)./(pi.*(d.^3)));%Maxium Shear Stress (in Pa)

Taua=KB.*((8.*Fa.*D)./(pi.*(d.^3)));%Alternating Shear Stress Component (in Pa)

Taum=KB.*((8.*Fm.*D)./(pi.*(d.^3)));%Midrange Shear Stress Component (in Pa)

Sty=[250
245
240
235
230
228
226];% Tensile Yield Strength (in kpsi)

Sty=Sty./0.000000145038; %Tensile Yield Strength (in Pa)

Sut=[275
270
265
260
255
253
251];% Ultimate Tensile Strength (in kpsi)

Sut=Sut./0.000000145038;% Ultimate Tensile Strength (in Pa)

Ssy=0.45.*Sty;%Shear Yield Strength (in Pa)

Sse=Ssa.*(1./(1-(Ssm./Sut)));%Shear Endurance Strength [Modified Goodman Line] (in Pa)

nmax=Ssy./Taumax;%Factor of safety  against Shear Yield; with respect to maximum Shear Stress

nf=sqrt(1./((Taua./Sse).^2+(Taum./Ssy).^2))%Factor of safety against fatigue failure [ASME Ellipitic Criterion]

Na=(d.^4.*G)./(8*k.*(D.^3));%Number of Active Coils

Nt=Na+2;%Total number of Coils for Squared and Grounded Ends

Ls=d.*Nt;%Solid Length (in m)

Fs=(1+Xi)*Fmax;%Force at Solid Length (in N)

Taus=KB.*((8*Fs.*D)./(pi.*(d.^3)));%Shear Stress at Solid Length (in Pa)

ns=Ssy./Taus;%Factor of safety against Shear Yield; with respect to Shear Stress at solid length 

Lcr=(pi.*D).*(sqrt(2*(E-G)./(2*G+E)));%Critical Length for Buckling; both ends pivoted

W=1/4*gamma*pi^2.*d.^2.*Nt.*D;%Weight of Spring (in N)

fom=-RC.*W;%Figure of Merit
