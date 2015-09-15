% demoltmcompare
close all

load ltmvary

fz = 14;
obj = 0;  %obj=1 slide, obj=0 paper


if obj == 1
 cm = jet(6);
 linestylelist= {'-','-','-','-','-','-'};
else
 cm = zeros(6,3);
 linestylelist= {'-','--','-.','d','x','s'};
end



figure
hold on

h(1)= plot(vary6)
h(2)= plot(vary10)
h(3)= plot(vary20)
h(4)= plot(vary30)

for i = 1:4

set(h(i),'color',cm(i,:) );
set(h(i),'linestyle',linestylelist{i});

end


axis([1 100 0 1])
set(h,'linewidth',2)
grid on
box on 
xlabel('iteration','fontsize',fz)
ylabel('normalized variance','fontsize',fz)
legend('3-cycle','5-cycle','10-cycle','15-cycle')
title('Comparison of LTM cycles','fontsize',fz);
set(gca,'fontsize',fz);
