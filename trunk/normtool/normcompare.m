%
% This routine compares the weighting of different norms
%
%
%
%
n = 10;

% s for Sobolev Norm 
  slist = [0 -0.5 -1]
% p for Sobolev Norm
  p     = 2;

% diffuse norm k
  k  = 1e-2;

% mix-norm weighting (pre-calculated using lkfind.m)
  load lk2000




[wsob ] = sobolevnorm(n,'f',p,slist);
[wdiff] = diffusenorm(n,'f',k);

w = [ lk(1:n+1);
      wdiff;
      wsob     ];

h = plot(0:n,w')
set(h,'linewidth',2);
grid on
legend('mix-norm','diffuse-norm','H_0','H_{-1/2}','H_{-1}')
xlabel('wave number');
ylabel('weight');


