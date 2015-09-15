%
% This routine compares the weighting of different norms
%
%
%
%
fz = 14;
fzl = 16;
obj = 0;  %obj=1 slide, obj=0 paper


if obj == 1
 cm = jet(6);
 linestylelist= {'-','-','-','-','-','-'};
else
 cm = zeros(6,3);
 linestylelist= {'-','--','-.','d','x','s'};
end


n = 10;

% s for Sobolev Norm 
  slist = [0 -0.5 -1]
% p for Sobolev Norm
  p     = 2;

% diffuse norm k
  k  = 0.01;

% mix-norm weighting (pre-calculated using lkfind.m)
  load lk2000




[wsob ] = sobolevnorm(n,'f',p,slist);
[wdiff] = diffusenorm(n,'f',k);

w = [ lk(1:n+1);
      wdiff;
      wsob     ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
figure 
hold on
box on
grid on
w = w([1,2,4,5,3],:);

for i = 1:5
    h(i) = plot(0:n,w(i,:));
    set(h(i),'linewidth',2,'color',cm(i,:),'linestyle',linestylelist{i},'markersize',8);
 
end




grid on
legend(h,'mix-norm','diffusion weight','H_{-0.5}','H_{-1}','L_{2}')
xlabel('wave number','fontsize',fzl,'interpreter','latex');
ylabel('weight','fontsize',fzl,'interpreter','latex');

set(gca,'fontsize',fz)



for i = 4:5
 gh(i) = plot(0:n,w(i,:),'linewidth',1,'color',cm(i,:));
end













