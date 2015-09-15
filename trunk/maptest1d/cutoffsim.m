function  [dist]= cutoffsim(S)


close all
figure 
hold on

colorlist = 'gcmybgcmybgcmyb'


munlist(1) = 1;
for i = 2:6
  temp = feval(S,munlist(i-1),'i');
  munlist(i) = temp(1);
end

y = 0:1/20000:1;
w0 = feval(S,y,'u');
p = 3

dist = zeros(length(munlist),length(munlist)+p);

for i=1:length(munlist)
  mun = munlist(i);
  
  w = 0*y;
  k = 0;
  %while k<7%w((length(y)+1)/2)==0;  
  while k< i+p  
    w = wfunS(y,k,mun,S);
    wd = w-w0;
    wd(1)=wd(2);, wd(end)=wd(end-1);
    dist(i,k+1) = trapz(y,abs(wd)/2);
    dist(find(dist>1))=1;
    h= stairs(y,w);
    % h = stairs(2/pi*asin(sqrt(y)),pi*sin(pi*y).*w);
    % h = stairs(sin(pi*y/2).^2,w);
    set(h,'color',colorlist(i));
    set(h,'linewidth',2); 
    k = k+1;

  end
end

h= plot(y, w0,'r');
set(h,'linewidth',2,'linestyle','--')
axis([0 1 0 10 ])
box on
grid on
xlabel('x')
ylabel('\omega(x)')
%h = plot(y,1./2./y);
%set(h,'color','b','linewidth',2,'linestyle','--')


figure
%h=  plot(0:size(dist,1)-1,dist')
h=  plot(0:size(dist,2)-1,dist')
set(h,'linewidth',2,'color','b')
box on
grid on
xlabel('k')
ylabel('\nu_{\mu_n}^k')

mapname = func2str(S)
mapname = mapname(1:end-3);

title(sprintf('%s map cutoff',mapname))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
h=  semilogy(0:size(dist,2)-1,dist')
set(h,'linewidth',2,'color','b')
box on
grid on

xlabel('k')
ylabel('log\nu_{\mu_n}^k') 
title(sprintf('%s map cutoff in log scale',mapname))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cutoffnormalize(dist,0.5)
title(sprintf('%s map cutoff normalized',mapname))
xlabel('normalized k')
ylabel('log\nu_{\mu_n}^k') 
grid on
box on



