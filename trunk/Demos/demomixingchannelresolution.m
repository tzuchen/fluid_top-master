

load Sdata
load Mexact

fz  = 14;
dx  = 0.02;

% data set 1
  load M200Markov
  load Mexact
  Ma = Mm;  

% data set 2
  load Sdata
  %load Mexact_high
 
 


close all
%Ma=Mm;
%Mb = Mc;
Mb(150:200) = Mb(165:215);

Mf{4} = Mb;
Mf{1} = M100;
Mf{2} = M200;
Mf{3} = M400;

w=1;
framelist = [10 50 100 200];
pl        = [1 2 3 4];
 for k = 1:4
    for i = 1:length(framelist)
 
   Maf = frame2im(Mf{k}(framelist(i)));
   
  

     [l,b]=ind2sub([4,4],i);
     axes('position',[(i-1)*0.25+dx ((5-k-1)*0.25)+dx  0.25-2*dx  0.25-2*dx])
 % j
 %  if(j==2)
 %     Maf;
 %  end

    if(k==4) 
       imagesc(Maf(:,:,3));
    else
       imagesc(Maf);
    end
    axis equal
    axis tight
    set(gca,'XTick',[],'XTickLabel',[]);
    set(gca,'YTick',[],'YTickLabel',[]);
    colormap gray
    mytitle = sprintf('(%s) k = %d',char(96+w),framelist(i));
    text(0, -12, mytitle,'fontsize',fz);
w=w+1;
   %  [l,b]=ind2sub([4,4],8+i);
   %  axes('position',[(i-1)*0.25+dx  ((k-1)*0.25)+dx  0.25-2*dx  0.25-2*dx]);

    end
end
