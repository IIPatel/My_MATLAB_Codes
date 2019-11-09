clear all; close all; clc;

x=[100:50:400]';
y=[10.63 13.03 15.04 16.81 18.42 19.90 21.27]';

xx=218;

p=(xx-x(1))/50;

del_n_y=zeros(7,7);
del_n_y(:,1)=y;
i=1;
for c=2:length(y)
    for r=1:(7-i)
        del_n_y(r,c)=del_n_y(r+1,c-1)-del_n_y(r,c-1);
    end
    i=i+1;
end

coeff1=1;
interpol=y(1);

for i=2:6
    coeff1=(p-(i-2))*coeff1;
    denom=factorial(i-1);
    interpol=interpol+del_n_y(1,i)*coeff1/denom;
end

interpol