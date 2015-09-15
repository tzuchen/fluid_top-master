% This is the demo of why cutoff happen
%
% Tzu-Chen Liang  9-27-2007

close all
figure
hold on

x = -3:0.001:3;
mu = 0.5;
sigma = 0.3;

y = normpdf(x,mu,sigma);
h = plot(x,y);
set(h,'linewidth',2);

y = normpdf(x,mu-2,sigma);
h = plot(x,y);
set(h,'linewidth',2);

y = normpdf(x,mu-0.3,sigma);
h = plot(x,y,'--');
set(h,'linewidth',2);
box on

