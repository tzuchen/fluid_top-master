%test logistictent.m
%
% This is the routine to check logistic map is
% actually equivalent to tantmap
%
% Tzu-Chen Liang 5-17-2007


y = 0:1/10000:1;

close all
figure 
hold on
colorlist = 'bgcmybgcmybgcmy'


clear munlist
%nlist =n0*[ ];
munlist(1) = 1;
for i = 2:7
  munlist(i) = munlist(i-1)/2
end

y = 0:1/10000:1;
w0 = 0*y+1;
p = 3

dist = zeros(length(munlist),length(munlist)+p);

for i=1:length(munlist)
  mun = munlist(i);
  
  w = 0*y;
  k = 0;
    
  while k< i+p  
    w = wfuntent(y,k,mun);

    wd = w-w0;
    dist(i,k+1) = trapz(y,abs(wd)/2);

    wlog = 1./pi./sqrt(y.*(1-y)).*w;
    %stairs(y,wlog);
    h= stairs(y,w);
    set(h,'color',colorlist(i));
    set(h,'linewidth',2); 
    k = k+1;

  end
end


h= plot(y, w0,'r');
set(h,'linewidth',2,'linestyle','--')
axis([0 1 0 10 ])
box on

h = plot(y,1./2./y);
set(h,'color','b','linewidth',2,'linestyle','--')


figure
h=  plot(dist')
set(h,'linewidth',2)
figure
h=  semilogy(dist')
set(h,'linewidth',2)



