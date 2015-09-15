function dist= testfunubound(n0)
%
%
%
close all
figure 
hold on

colorlist = 'gcmybgcmybgcmyb'



%nlist =n0*[ ];
munlist(1) = 1;
for i = 2:8
  munlist(i) = (1-sqrt(1-munlist(i-1)))/2
end

y = 0:1/20000:1;
w0 = 1./(pi*(y.*(1-y)).^0.5);
p = 3

dist = zeros(length(munlist),length(munlist)+p);

for i=1:length(munlist)
  mun = munlist(i);
  
  w = 0*y;
  k = 0;
  %while k<7%w((length(y)+1)/2)==0;  
  while k< i+p  
    w = wfun(y,k,mun);
    wd = w-w0;
    wd(1)=wd(2);, wd(end)=wd(end-1);
    dist(i,k+1) = trapz(y,abs(wd)/2);
    dist(find(dist>1))=1;
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
title('logistic map cutoff')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
h=  semilogy(0:size(dist,2)-1,dist')
set(h,'linewidth',2,'color','b')
box on
grid on

xlabel('k')
ylabel('log\nu_{\mu_n}^k') 
title('logistic map cutoff in log scale')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function y = S(x)
%y = 4*x.*(1-x);
%function y = Si(x)
%y = (1-sqrt(1-x))/2;
%function y = Sp(x)
%y = 8-4*x;
%function y = f(x)
%y = (1-sqrt(1-4*x.^2/pi^2))/2;
%function y =wn(x);
%y = 1./(pi*sqrt(x.*(1-x)));







