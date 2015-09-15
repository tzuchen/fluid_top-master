% n-D^star chart
% 
% This routine produce the relationship of n(number of grids) 
% and D^*(effective D) chart, and shows they are linear
%
% Tzu-Chen Liang 7/28/2006

clear D mapfunction param
nlist = [50:100:1500];

mapfunction{1} = @standardmap;
param{1}{1}    = 0.8;

mapfunction{2} = @standardmap;
param{2}{1}    = 0.3;

mapfunction{3} = @twistmap;
param{3}{1}    = 0.5;
param{3}{2}    = 0.5;
param{3}{3}    = 1;

mapfunction{4} = @bakersmap;
param{4}{1}    = [];

mapfunction{5} = @arnoldscatmap
param{5}{1}    = 0;

for caseno = 1:size(mapfunction,2);

   for i = 1:length(nlist)
       n      = nlist(i);
       D(caseno,i) = Destimate(mapfunction{caseno},param{caseno},n)

   end
end

 save /home/tzuchen/Desktop/T'emp computation result'/nddata nlist D
 load /home/tzuchen/Desktop/T'emp computation result'/nddata nlist D

close all
figure 
hold on
h = loglog(1./(nlist.^2),D')
set(h,'linewidth',2)
h = loglog(1./(nlist.^2),D','x')
set(h,'linewidth',2)
set(gca,'YScale','log')
set(gca,'XScale','log')
legend('standard map, 0.8','standard map, 0.3','twistmap','bakers map','arnold cap map','Location','NorthWest')
xlabel('1/n^2')
ylabel('D^*')
grid on
box on

