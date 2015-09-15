% Standard Map cutoff plot 
% epsilon = 0.3 (I think)
% initial function  cos(2pi x) i think
%
% Tzu-Chen Liang 10-17-2007

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

load mixnorm_2500_60s
load mixnorm_5000_60s
load mixnorm_10000_60s
load mixnorm_20000_60s
load mixnorm_40000_60s
load mixnorm_80000_60s


nlist = 2*[mixnorm_2500_60s(:,2),mixnorm_5000_60s(:,2),mixnorm_10000_60s(:,2),mixnorm_20000_60s(:,2),mixnorm_40000_60s(:,2),mixnorm_80000_60s(:,2)]';
[dl,iter] = size(nlist);

dlist = [2500 5000 10000 20000 40000 80000];

for i = 1:dl
  w{i} = sprintf('$n = %d^2$',dlist(i));
end



close all
figure
hold on

for i = 4:6
thline = plot(0:50,nlist(i,1:51),'color',[0 0 0],'linewidth',1);
end

for i = 1:6

h(i) = plot(0:50,nlist(i,1:51));
set(h(i),'color',cm(i,:) );
set(h(i),'linestyle',linestylelist{i});
end
 

axis([0 40 0.4 1.015]);
set(h,'linewidth',2);
xlabel('$k$ (iteration)','fontsize',fzl,'interpreter','latex');
ylabel('var$(f)$','fontsize',fzl,'interpreter','latex');
box on
grid on
legend(h,w,'interpreter','latex');
title('Simulation of Standard Map with extra smoothing steps','fontsize',fzl,'interpreter','latex');
set(gca,'fontsize',fz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%noval = 0.7
noval = (1+nlist(end,end))/2

for i = 1:dl
  for j = 1:iter-1
    if nlist(i,j+1)>=nlist(i,j);
      nlist(i,j+1) = nlist(i,j)-1e-9; 
    end
  end
end


for i = 1:dl
  itern(i)   = interp1(nlist(i,:),0:iter-1,noval);
end



figure
hold on

if obj==1
plot([0 1],[1 1],'linewidth',2,'linestyle','--');
plot([1 1],[1 nlist(end:end)],'linewidth',2,'linestyle','--');
plot([1 3],[nlist(end:end) nlist(end:end)],'linewidth',2,'linestyle','--');
end

for i = 1:dl

h(i) = plot([0:iter-1]/itern(i),nlist(i,:));
set(h(i),'color',cm(i,:));  
set(h(i),'linestyle',linestylelist{i});
end
legend(h,w,'interpreter','latex');
axis([0 3 0.4 1.015]);
set(h,'linewidth',2);

box on
grid on
xlabel('$k$/$t_n$ (normalized iteration)','fontsize',fzl,'interpreter','latex');
ylabel('var$(f)$','fontsize',fzl,'interpreter','latex');
title('Simulation of Standard Map with extra smoothing steps','fontsize',fzl,'interpreter','latex');
set(gca,'fontsize',fz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zoom 
vs= axis
delta=  0.1
%boxpos1 = [1.41,0.46];
%boxpos2 = [0.3,0.93];
boxpos1 = [1.47,0.48];
boxpos2 = [0.34,0.94];


aratio = [vs(2)-vs(1) vs(4)-vs(3)];


box1pt = [boxpos1(1) boxpos1(1)+delta*aratio(1) boxpos1(2) boxpos1(2)+delta*aratio(2)];
v1(1) = plot([box1pt(1) box1pt(1)],[box1pt(3) box1pt(4)]);
v1(2) = plot([box1pt(2) box1pt(2)],[box1pt(3) box1pt(4)]);
v1(3) = plot([box1pt(1) box1pt(2)],[box1pt(3) box1pt(3)]);
v1(4) = plot([box1pt(1) box1pt(2)],[box1pt(4) box1pt(4)]);
set(v1,'linewidth',2)


box2pt = [boxpos2(1) boxpos2(1)+delta*aratio(1) boxpos2(2) boxpos2(2)+delta*aratio(2)];
v2(1) = plot([box2pt(1) box2pt(1)],[box2pt(3) box2pt(4)]);
v2(2) = plot([box2pt(2) box2pt(2)],[box2pt(3) box2pt(4)]);
v2(3) = plot([box2pt(1) box2pt(2)],[box2pt(3) box2pt(3)]);
v2(4) = plot([box2pt(1) box2pt(2)],[box2pt(4) box2pt(4)]);
set(v2,'linewidth',2)



delta=1

axes('position',[0.62 0.31 0.2 0.2])
hold on
box on


for i = 1:dl
   h(i) = plot([0:iter-1]/itern(i),nlist(i,:));
   set(h(i),'color',cm(i,:));  
   set(h(i),'linestyle',linestylelist{i});
   if obj ==1
      set(h(i),'linewidth',2)
   else
      set(h(i),'linewidth',1)
   end
end
axis(box1pt)




delta=1

axes('position',[0.44 0.66 0.2 0.2])
hold on
box on


for i = 1:dl
   h(i) = plot([0:iter-1]/itern(i),nlist(i,:));
   set(h(i),'color',cm(i,:));  
   set(h(i),'linestyle',linestylelist{i});
   if obj ==1
      set(h(i),'linewidth',2)
   else
      set(h(i),'linewidth',1)
   end
end
axis(box2pt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear xi yi areai
%figure 
%hold on
%box on
xi = 0:1e-4:3;
for i = 1:dl

 yi(i,:) = interp1([0:iter-1]/itern(i),nlist(i,:),xi,'spline');
 areai(i) = trapz(xi,abs(yi(i,:)-noval));
end
 itern2 =itern;
 noval2 =noval;
 nlist2 =nlist;
 totalarea2 = abs(nlist(end,end)-noval)*3; 
 areai2 = areai;
 dlist2 = dlist;

%g1 = plot(log(log(dlist.^2)),1.6563-areai)
%set(g1,'linewidth',2);
%g2 = plot(log(log(dlist.^2)),1.6563-areai,'s')
%set(g2,'markersize',10,'linewidth',2);

