% demostandardmapfreqevolve

jlist = 1:15


for j=jlist
  fname = ['wnlog_40000_', num2str(j),'o']
  load(fname)
end

close all
figure 
hold on
box on

for j=jlist

    fname = ['wnlog_40000_', num2str(j),'o'];
    w     = eval(fname);
    h     = loglog(w(:,1),w(:,2));
    set(h,'color',[0 0 0])

end
set(gca,'xscale','log')
set(gca,'yscale','log')


xlabel('wave number','fontsize',16,'interpreter','latex');
ylabel('magnitude','fontsize',16,'interpreter','latex');
set(gca,'fontsize',fz)
%title('Standard Map, $\epsilon=0.3$','interpreter','latex')

