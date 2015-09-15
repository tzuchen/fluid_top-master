
x = 0.001:0.05:0.5; 
y = 0.5*ones(size(x));

param{1} = 0.5;
param{2} = 0.5;
param{3} = 0.1;
close all
figure
hold on

for i = 1:1000
[x,y] = twistmap(x,y,1,1,param);
plot(x,y,'.')
axis([0 1 0 1])

end
