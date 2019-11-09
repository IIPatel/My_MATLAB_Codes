clc
clear all
close all

%%Calling jpg file to trace geometry and obtain the coordinates%%

%  I = imread('images.jpg');
%  [BW, xi, yi] = roipoly(I); %returns the x- and y-coordinates of the polygon vertices in xi and yi.

%Copied coordinates from the above procdure%
x=[20.6848
   20.6848
   20.6848
   20.6848
   20.6848
   20.6848
   23.0735
   24.6659
   26.6564
   27.6517
   29.6422
   30.8365
   32.2299
   31.0355
   29.6422
   29.0450
   27.6517
   27.0545
   25.2630
   26.4573
   27.8507
   28.8460
   28.0498
   26.6564
   26.2583
   24.8649
   23.8697
   23.6706
   26.4573
   27.8507
   26.6564
   25.4621
   24.4668
   23.6706
   22.8744
   21.8791
   20.8839
   20.0877
   18.8934
   18.2962
   17.6991
   16.5047
   15.9076
   15.1114
   13.9171
   13.5190
   13.1209
   11.7275
   10.9313
    9.9360
    8.7417
    7.3483
    7.1493
    8.5427
    9.7370
   11.1303
   10.9313
    9.7370
    8.7417
    7.3483
    8.9408
    8.3436
    6.7512
    4.9597
    3.9645
    2.1730
    4.9597
    5.7559
    7.9455
   10.5332
   11.9265
   13.9171
   13.9171
   13.9171
   13.9171
   13.9171
   13.9171
   14.7133
   15.9076
   17.3009
   18.8934
   20.4858
   20.4858
   20.4858
   20.6848];

y=[-40.8081
  -38.8175
  -37.4242
  -36.2299
  -34.0403
  -33.0450
  -33.8412
  -33.8412
  -34.4384
  -34.4384
  -34.4384
  -34.4384
  -34.4384
  -32.2488
  -30.8555
  -29.8602
  -28.4668
  -27.2725
  -25.8791
  -25.2820
  -26.2773
  -25.8791
  -23.4905
  -22.0972
  -20.7038
  -19.1114
  -17.7180
  -16.7227
  -17.3199
  -16.9218
  -14.5332
  -12.5427
  -11.3483
  -10.3531
   -8.7607
   -7.3673
   -6.1730
   -5.1777
   -3.7844
   -2.3910
   -1.7938
   -1.7938
   -2.7891
   -3.5853
   -4.9787
   -6.5711
   -7.3673
   -9.1588
  -10.1540
  -11.5474
  -13.3389
  -14.9313
  -16.7227
  -15.5284
  -15.9265
  -15.9265
  -18.1161
  -19.9076
  -21.5000
  -24.2867
  -25.6801
  -27.2725
  -29.8602
  -31.4526
  -32.8460
  -35.0355
  -35.0355
  -33.8412
  -33.6422
  -33.2441
  -32.8460
  -32.8460
  -34.8365
  -36.0308
  -38.2204
  -39.0166
  -39.6137
  -41.4052
  -41.2062
  -41.2062
  -41.2062
  -41.2062
  -41.2062
  -40.8081
  -40.8081];

%Defining Pseudo Time Span
t=[0:1:length(x)-1]';
%Plots
subplot(231)
plot(x,y);%y=-y to correct the inverted image obtained by default
title('Alpine Tree Tracing')
xlabel('x-coordinates')
ylabel('y-coordinates')
subplot(232)
plot(t,x);
title('x against Pseudo time')
xlabel('Pseudo Time')
ylabel('x-coordinates')
subplot(233)
plot(t,y);
title('y against Pseudo time')
xlabel('Pseudo Time')
ylabel('y-coordinates')

% %X-Fourier Coefficients, for n=7 from Curve Fitting
%        a0 =       154.3;
%        a1 =     0.01341;
%        b1 =       11.35;
%        a2 =     0.07254;
%        b2 =      0.1954;
%        a3 =     -0.1179;
%        b3 =       1.145;
%        a4 =    -0.05841;
%        b4 =       2.357;
%        a5 =     0.07658;
%        b5 =     -0.1326;
%        a6 =    -0.06081;
%        b6 =       1.513;
%        a7 =    -0.2317;
%        b7 =       3.819;
%        w =      0.4189;
%        x=t;
%   ft= a0 + a1*cos(x*w) + b1*sin(x*w) + a2*cos(2*x*w) + b2*sin(2*x*w) + a3*cos(3*x*w) + b3*sin(3*x*w) + a4*cos(4*x*w) + b4*sin(4*x*w) + a5*cos(5*x*w) + b5*sin(5*x*w) + a6*cos(6*x*w) + b6*sin(6*x*w) + a7*cos(7*x*w) + b7*sin(7*x*w);
% X=ft;
% %Y-Fourier Coefficients, for n=7 from Curve Fitting
%        a0 =       146.7;
%        a1 =      -14.09;
%        b1 =     -0.5179;
%        a2 =      -1.629;
%        b2 =         0.1;
%        a3 =      -1.729;
%        b3 =      0.1386;
%        a4 =      -0.925;
%        b4 =    -0.01656;
%        a5 =      -2.068;
%        b5 =     -0.1326;
%        a6 =      0.5545;
%        b6 =     -0.1157;
%        a7 =       -4.77;
%        b7 =      0.0598;
%        w =      0.4189;
%  ft= a0 + a1*cos(x*w) + b1*sin(x*w) + a2*cos(2*x*w) + b2*sin(2*x*w) + a3*cos(3*x*w) + b3*sin(3*x*w) + a4*cos(4*x*w) + b4*sin(4*x*w) + a5*cos(5*x*w) + b5*sin(5*x*w) + a6*cos(6*x*w) + b6*sin(6*x*w) + a7*cos(7*x*w) + b7*sin(7*x*w);
% Y=-ft;
% subplot(234)
% plot(X,Y);
% title('Tree, n=7')
% xlabel('X-coordinates')
% ylabel('Y-coordinates')
% subplot(235)
% plot(t,X);
% title('X against Pseudo time')
% xlabel('Pseudo Time')
% ylabel('X-coordinates')
% subplot(236)
% plot(t,Y);
% title('Y against Pseudo time')
% xlabel('Pseudo Time')
% ylabel('Y-coordinates')
% 
%                
% 
% 
