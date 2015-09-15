
%[Astruc] = maprefine2(500,[],@standardmap,0.3)
%[Xa1,var2,var1,mnorm,snorm,sv,Mm,alln] = testAf(Astruc,[],'cosy',30,1,0,0)


flist= [5 1 6 2 10 3 16 4];
f=1;
for p = 0:3

for q = 0:1


   h = axes('Position',[ p/4+0*d 0.95*(q/2+0*d)  1/4-d 0.95*(1/2-d) ]);
  
  myf=Mm(flist(f)).cdata;
  imagesc(myf);
  colormap(gray)
  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  set(gca,'TickLength',[0 0]);
  grid off
s = sprintf('k=%d',flist(f)-1);
text(0,-15,s,'fontsize',fz)
  axis equal
  axis tight
  f=f+1
 
% if q==3
%  s = sprintf('Stroock''s Mixer');
%  %text(0,245,s,'fontsize',fz)
%  text(0,-45,s,'fontsize',fz)
% end

end

end




