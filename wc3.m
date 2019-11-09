clc
clear all
close all

o=250e-3/2;%Outlet half-width
i=o*2;%Inlet half-width;contraction ratio is i/o
Lc=600e-3;%Length of contraction section
xc=[0:1e-3:Lc]';%x coordinate points of section
yc=i-(i-o).*(-2.*(xc./Lc).^3+3.*(xc./Lc).^2);%y coorinates to model 3rd order curve
zc=zeros(601,1);%z coordinates to model 5th ordercurve
subplot(223)
plot(xc,yc)
title('Contraction Section 3rd order (2:1)')
hold on

Lt=0.8;%Length of contraction section
xt=[0.6:1e-3:(Lc+Lt)]';%x coordinate points of section
for i=1:length(xt)
    yt(i)=o;%y coordinates
    i=i+1e-3;
end
yt=yt';
zt=zeros(801,1);%z coordinates to model 5th ordercurve

%  filename = '3poly Contraction Section.xlsx';
%  sheet = 1;
%  xlRange = 'A2';
%  xlswrite(filename,xc,sheet,xlRange)
%  filename = '3poly Contraction Section.xlsx';
%  sheet = 1;
%  xlRange = 'B2';
%  xlswrite(filename,yc,sheet,xlRange)
%   filename = '3poly Contraction Section.xlsx';
%  sheet = 1;
%  xlRange = 'C2';
%  xlswrite(filename,zc,sheet,xlRange)
% 
%  