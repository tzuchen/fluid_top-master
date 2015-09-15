close all
figure
hold on
box on

n= 60
x(1) = 0.87;
for i=1:n
  x(i+1)=x(i)*0.93;
end


x = x+ 0.1*rand(n+1,1);
x = sort(x);
x = x(end:-1:1);
x = x/2+0.5;

h = stairs([0:n],x)
set(h,'linewidth',2)

axis([0 n+10 0.5 1])

M = sum(x-0.5)
Mn = ceil(M/(max(x)-0.5));
xM = 0.5*ones(n+1,1);
xM(1:Mn)=max(x);

g=  stairs(xM)
set(g,'linewidth',2,'color','r')
