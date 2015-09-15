function  [dist]= funcutoffsim(S)


dx = 1/10000;
x = dx/2:dx:1-dx/2;
xs = x;

close all
figure
%hold on
k = 0
for i = [5 6 7 8 1 2 3 4]

[l,b]=ind2sub([4,2],i)
axes('position',[(l-1)*0.25+0.02  ((b-1)*0.5)+0.02  0.20  0.40])

 
 h = plot(xs,sin(2*pi*x));
 set(h,'linewidth',2)
 set(gca,'XTick',[],'XTickLabel',[])
 set(gca,'YTick',[],'YTickLabel',[])
 axis([0 1 -1 1])
 box on
 title(sprintf('k = %d',k))
 x = feval(S,x);
 k = k +1;
end


dist= [];

