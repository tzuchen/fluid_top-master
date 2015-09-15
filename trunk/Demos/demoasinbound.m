%demoasinbound.m

close all
figure 
hold on
nlist =[1e5 1e10 1e15 1e20];
w = 4;


for i=1:length(nlist)
    n = nlist(i);
    k = 1:log(n)/log(w)+1;
    p = w.^k/n;
    f = 1-4/pi*asin(p.^0.5);
    %f(find(f<0)) =0;
         

    h = plot(real(f))
    set(h,'linewidth',2)
end

for i = 1:length(nlist)
  wl{i} = sprintf('n = %g',nlist(i));
end

box on
grid on
axis([0 40 0 1])
title('lower bounds as a function of n')
xlabel('k')
ylabel('lower bounds')
legend(wl)
