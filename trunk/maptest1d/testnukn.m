%testnukn

close all
figure 
hold on

k = 0:15
nlist = 1:2:15

for i=1:length(nlist)
   n =nlist(i);
   h= plot(k,max(1-2.^(1+k-n),0));
   set(h,'linewidth',2)
end

grid on
box on
axis([0 15 0 1])
xlabel('k')
ylabel('\nu_{\mu_n}^k')
title('Cutoff of tent map')
