%demo of standardmap cutoff
%now the initial condition is a point

load mixnorm_2500_100ic1
load mixnorm_5000_100ic1
load mixnorm_10000_100ic1
load mixnorm_20000_100ic1
load mixnorm_40000_100ic1
load mixnorm_80000_100ic1

mixnorm_2500_50 = mixnorm_2500_100ic1(1:50,:);
mixnorm_5000_50 = mixnorm_5000_100ic1(1:50,:);
mixnorm_10000_50 = mixnorm_10000_100ic1(1:50,:);
mixnorm_20000_50 = mixnorm_20000_100ic1(1:50,:);
mixnorm_40000_50 = mixnorm_40000_100ic1(1:50,:);
mixnorm_80000_50 = mixnorm_80000_100ic1(1:50,:);

p=2;
nlist = 2*[mixnorm_2500_50(:,p),mixnorm_5000_50(:,p),mixnorm_10000_50(:,p),mixnorm_20000_50(:,p),mixnorm_40000_50(:,p),mixnorm_80000_50(:,p)]';
[dl,iter] = size(nlist);

dlist = [2500 5000 10000 20000 40000 80000];

for i = 1:dl
  w{i} = sprintf('n = %d',dlist(i));
end



close all
figure
hold on
h = plot(nlist')
set(h,'linewidth',2);
xlabel('iteration')
ylabel('var(f)')
box on
grid on
legend(w);
title('Simulation of Standard Map with different number of grids')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noval = 0.7

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

clist = 'bgrcmyk'

figure
hold on

for i = 1:dl

h(i) = semilogy([0:iter-1]/itern(i),nlist(i,:));
set(h(i),'color',clist(i));  
end
legend(w);
%axis([0 5 0.43 1])
set(h,'linewidth',2);

box on
grid on
xlabel('nornalized iteration')
ylabel('total variation distance')
title('Simulation of Standard Map with different number of grids')
