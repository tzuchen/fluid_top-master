function testbound(n)
%testbound.m
%
% To test two lower bounds
% for pi measure and uniform measure respectively
% of the convergence trajectory.
%
%
% Tzu-Chen Liang 4-19-2007


 dx   = 1/n;
 x    = dx/2:dx:1-dx/2;  
 xinv = (1/pi./(x.*(1-x)).^0.5)';
 xinv = xinv/sum(xinv);


A = maprefine1d2(n,@logisticmap);
%[x,mixnorm1,mixnorm2]=Atest1d(A,'cosx',30);
[x,m1,m2,m3]=Atest1d(A','sing',30,xinv);

 w = 4 % for logistic map
 k = 1:log(n)/log(4);

 p = w.^k/n;
 f = 1-4/pi*asin(p.^0.5);


k = 1:100;
lb = 1-(w.^(k))/n;
lb2 = 2.^(-k)*n-1;


close all
figure
hold on

plot(m1);
%plot(m2);
plot(m3);
plot(f,'y');
plot(lb,'r');

%plot(lb2,'g');



axis([1 30 0 1.1]);
