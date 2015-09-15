% This routine generates the figures for Standardmap example
%
%
% 
% Tzu-Chen Liang   2-19-2007

%load Mm

close all

Ma= Mm;

framelist = [40 60 80 100 1 2 5 20 ]

for i = 1:length(framelist)
   Maf = frame2im(Ma(framelist(i)));
  
    [l,b]=ind2sub([4,2],i)
    axes('position',[(l-1)*0.25+0.02  ((b-1)*0.25)+0.02  0.20  0.20])
    imagesc(Maf)
    axis equal
    axis tight
    set(gca,'XTick',[],'XTickLabel',[])
    set(gca,'YTick',[],'YTickLabel',[])
    colormap gray
     title(sprintf('t = %d',framelist(i)))
 
end
