% This routine generates the figures for Mixing channel example
%
%
% 
% Tzu-Chen Liang    2-4-2007
%                 10-17-2007

fz  = 14;
dx  = 0.02;

% data set 1
  load M200Markov
  load Mexact
  Ma = Mm;  

% data set 2
  load Sdata
  %load Mexact_high
  Ma(100:300)=M200(95:295);
  %Mb = Mc;


close all
%Ma=Mm;
%Mb = Mc;



framelist = [1 51 151 251];
pl        = [1 2 3 4];

for i = 1:length(framelist)
   Maf = frame2im(Ma(framelist(i)));
   Mbf = frame2im(Mb(framelist(i)));
  

     [l,b]=ind2sub([4,2],i);
     %axes('position',[(l-1)*0.25+dx ((b-1)*0.25)+dx  0.25-2*dx  0.25-2*dx])
     axes('position',[(l-1)*0.25+dx ((b-1)*0.5)+dx  0.25-2*dx  0.5-2*dx])

    imagesc(Maf);
    axis equal
    axis tight
    set(gca,'XTick',[],'XTickLabel',[]);
    set(gca,'YTick',[],'YTickLabel',[]);
    colormap gray
    mytitle = sprintf('(%s) k = %d',char(96+4+pl(i)),framelist(i)-1);
    text(0, -18, mytitle,'fontsize',fz);

     [l,b]=ind2sub([4,2],4+i);
     %axes('position',[(l-1)*0.25+dx  ((b-1)*0.25)+dx  0.25-2*dx  0.25-2*dx]);
     axes('position',[(l-1)*0.25+dx  ((b-1)*0.5)+dx  0.25-2*dx  0.5-2*dx]);

    imagesc(Mbf(:,:,3));
    axis equal
    axis tight
    set(gca,'XTick',[],'XTickLabel',[]);
    set(gca,'YTick',[],'YTickLabel',[]);
    colormap gray
    mytitle = sprintf('(%s) k = %d',char(96+pl(i)),framelist(i)-1);
    text(0 ,-18,mytitle,'fontsize',fz);
    
end
