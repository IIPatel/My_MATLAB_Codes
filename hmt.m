clc
clear

x=[0:15:105]'

Case1=[58 53.5 53.5 49 38.8 32 31 29.4]';

Case2=[66.8 60.7 60.6 55 42 33.6 32.2 30.2]';

Case3=[80 71.7 71.5 64.3 47.2 35.9 34 31.3]';

Case4=[89.7 80.7 80.6 72.1 52 38.5 35.8 32.5]';

f=fittype('poly1')
F1=fit(x, Case1, f)
F2=fit(x, Case2, f)
F3=fit(x, Case3, f)
F4=fit(x, Case4, f)


plot(F1, x, Case1)
hold on
plot(F2, x, Case2)
hold on
plot(F3, x, Case3)
hold on
plot(F4, x, Case4)