%Demostandardmapevolve
%
%
% Tzu-Chen Liang 5-30-2007




%[Astruc] = maprefine2(1000,[],@standardmap,0.3)


for casen = 2:2

 if casen==1
   A = Astruc.A;
   ic = 'cosx';
 else
   A = Astruc.A'; 
   ic = 'sing';
 end

figure
hold on

klist = [1 3 6 11 21 31]

 for i = 1: length(klist)
  [Xa{i},var2,var1] = testAs(A,[],0,ic,klist(i),0);
 end


q = 0;
for i = [4 5 6 1 2 3]
 q = q+1;
 [l,b]=ind2sub([3,2],i)
 axes('position',[(l-1)/3+0.02  ((b-1)/2)+0.02  0.3  0.45])
 imagesc(Xa{q});

 box on          
 axis equal  
 axis tight       
 set(gca,'ydir','normal');
 set(gca,'XTick',[],'XTickLabel',[])
 set(gca,'YTick',[],'YTickLabel',[])
 title(sprintf('k = %d',klist(q)-1))
 hold on

end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 
