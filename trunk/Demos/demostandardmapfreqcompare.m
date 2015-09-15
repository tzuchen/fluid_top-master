%demostandardmapfreqcompare

close all
noffile = 20;
fz = 14;
fzl = 16;

jlist = [3 7 20]
obj = 0

if obj == 1
 cm = lines(6);
 linestylelist= {'-','-','-','-','-','-'};
else
 cm = zeros(6,3);
 linestylelist= {'-','-','-','d','x','s'};
end



for j=jlist
  

  fname = ['wnlog_40000_', num2str(j),'o'];
  load(fname)
  fname = ['wnlog_40000_', num2str(j),'s'];
  load(fname)
  fname = ['wnlog_40000_', num2str(j),'f'];
  load(fname)

end


figure 
hold on
box on
flist = floor([10^3.61 10^3.8 10^4]);

k = 1
for j = jlist

    fname = ['wnlog_40000_', num2str(j),'o'];
    w     = eval(fname);
    h(k) = loglog(w(:,1),w(:,2));
    if obj==0
      mk(1) = plot(w(flist(k),1),w(flist(k),2),'s','markersize',8,'linewidth',2,'color',cm(1,:)); 
    end

    fname = ['wnlog_40000_', num2str(j),'f'];
    w     = eval(fname);
    g(k) = loglog(w(:,1),w(:,2));
    if obj==0 
      mk(2) = plot(w(flist(k),1),w(flist(k),2),'d','markersize',8,'linewidth',2,'color',cm(2,:)); 
    end

    fname = ['wnlog_40000_', num2str(j),'s'];
    w     = eval(fname);
    f(k) = loglog(w(:,1),w(:,2));
    if obj==0 
      mk(3) = plot(w(flist(k),1),w(flist(k),2),'x','markersize',8,'linewidth',2,'color',cm(3,:)); 
    end
 
    set(h(k),'linestyle',linestylelist{1},'color',cm(1,:));
    set(g(k),'linestyle',linestylelist{2},'color',cm(2,:));
    set(f(k),'linestyle',linestylelist{3},'color',cm(3,:));


    k = k+1;    
end





set(gca,'xscale','log');
set(gca,'yscale','log');
axis([1 10^5 10^-6 2])
if obj==0
  legend(mk,'Numerical diffusion',...
            'Physical diffusion',...
            'Smoothing step')
else
  legend([h(1),g(1),f(1)],'Numerical diffusion',...
            'Physical diffusion',...
            'Smoothing step','interpreter','latex')
end


xlabel('wave number','fontsize',fzl,'interpreter','latex')
ylabel('magnitude','fontsize',fzl,'interpreter','latex')

set(gca,'fontsize',fz)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


