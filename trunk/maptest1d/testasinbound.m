%testasinbound
%
%
%
close all
figure 
hold on
r= [5 ];
w = 4;

lbplot = 1;
ubplot = 0;
simuplot = 1;
explb  =  0;

nlist = 10.^r;
nit = 30;
%simulation

if simuplot == 1;
  for i=1:length(r)
     n = nlist(i);
     A = maprefine1d2(n,@logisticmap);
     [x,m1,m2,m3]=Atest1d(A','sing',10*nit);
     [x,m1,m2,m3]=Atest1d(A','sing',nit,x);
     h= plot(m3,'g');
     set(h,'linewidth',2)
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lower bound
if lbplot == 1
  for i=1:length(r)
    n = nlist(i);
    k = 1:log(n)/log(w);
    p = w.^k/n;
    f = 1-4/pi*asin(p.^0.5);
    f(find(f<0)) =0;
    
    h = plot(f)
    set(h,'linewidth',2)
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear f
%upper bound

if ubplot ==1
  for i=1:length(r)
    n = nlist(i);
    xr = 1/n;
    clear f
    indend = 0;
    for j = 1:30
      p =1-xr;

      f(j) = 2/pi*asin(p.^0.5);
      if indend==0
         xr = logisticmap(xr);
      end
      if xr>0.5
         indend = 1;       
      end

   
    end
    plot(f,'r')
    plot(f,'rx')

  end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%exponential upper bound
if explb == 1;
  w = 0.24; 
  for i=1:length(r)
   n = nlist(i);
   k = 0.7383*log(n);

   kl = 0:k:60;
   myfac = 1-n^(0.7383*log(w)+1);
   p = zeros(1,length(kl))
   p(1) = 1;
   for j = 2:length(kl)
     p(j) = p(j-1)*myfac;
   end
   %plot(p)
   h= plot(kl,p)
   set(h,'color','y')
   plot(kl,p,'x')
  end
end

box on
grid on
xlabel('iteration')
ylabel('|x^k-\pi|/2')


legend('simulation','lower bound')
title('a lower bound, n=1e5')
axis([1 nit 0 1 ])
