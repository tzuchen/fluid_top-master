%testMapplot
%
%
% Tzu-Chen Liang 5-29 2007
close all



dx = 1/10000;
x = 0:dx:1;

figure
hold on
h1 = plot(x,tentmap(x));
set(h1,'linewidth',2)
grid on 
box on
%title('Tent Map')


%figure
%hold on
h1 = plot(x,logisticmap(x));
set(h1,'linewidth',2)
grid on 
box on
%title('Tent Map')

title('Tent Map and Logistic Map')


%figure
%hold on
h1 = plot(x,halfsinmap(x));
set(h1,'linewidth',2)



%figure
%hold on
h1 = plot(x,circlemap(x));
set(h1,'linewidth',2)

