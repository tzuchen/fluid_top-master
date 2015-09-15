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



framelist = [ 150 200 250 300 1 10 50 100];
%framelist = [ 100 140 180 200 1 10 30 70];
pl        = [5 6 7 8 1 2 3 4];

for i = 1:length(framelist)
   Maf = frame2im(Ma(framelist(i)));
   Mbf = frame2im(Mb(framelist(i)));
  

     [l,b]=ind2sub([4,4],i);
     axes('position',[(l-1)*0.25+dx ((b-1)*0.25)+dx  0.25-2*dx  0.25-2*dx])


    imagesc(Maf);
    axis equal
    axis tight
    set(gca,'XTick',[],'XTickLabel',[]);
    set(gca,'YTick',[],'YTickLabel',[]);
    colormap gray
    mytitle = sprintf('(%s) t = %d',char(96+8+pl(i)),framelist(i));
    text(0, -14, mytitle,'fontsize',fz);

     [l,b]=ind2sub([4,4],8+i);
     axes('position',[(l-1)*0.25+dx  ((b-1)*0.25)+dx  0.25-2*dx  0.25-2*dx]);

    imagesc(Mbf(:,:,3));
    axis equal
    axis tight
    set(gca,'XTick',[],'XTickLabel',[]);
    set(gca,'YTick',[],'YTickLabel',[]);
    colormap gray
    mytitle = sprintf('(%s) t = %d',char(96+pl(i)),framelist(i));
    text(0 ,-14,mytitle,'fontsize',fz);
    
end
