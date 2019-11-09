clc
clear all
close all

o=250e-3/2;%Outlet half-width
i=o*2;%Inlet half-width;contraction ratio is i/o
Lc=600e-3;%Length of contraction section
xc=[0:1e-3:Lc]';%x coordinate points of section
yc=i-(i-o).*(-20*(xc./Lc).^7+70.*(xc./Lc).^6-84.*(xc./Lc).^5+35.*(xc./Lc).^4);%y coorinates to model 3rd order curve
zc=zeros(601,1);%z coordinates to model 5th ordercurve
subplot(223)
plot(xc,yc)
title('Contraction Section 7th order (2:1)')

 filename = '7poly Contraction Section.xlsx';
 sheet = 1;
 xlRange = 'A2';
 xlswrite(filename,xc,sheet,xlRange)
 filename = '7poly Contraction Section.xlsx';
 sheet = 1;
 xlRange = 'B2';
 xlswrite(filename,yc,sheet,xlRange)
  filename = '7poly Contraction Section.xlsx';
 sheet = 1;
 xlRange = 'C2';
 xlswrite(filename,zc,sheet,xlRange)

 