function [X3, area] = polyintersect(X1,X2)
% given 2 sets of points, this function returns
% their intersection (as a polygon)
%
% Tzu-Chen Liang


[A1,b1]=pt2ab(X1(1,:),X1(2,:));
[A2,b2]=pt2ab(X2(1,:),X2(2,:));


A = [A1;A2];
b = [b1;b2];

ang  = [cart2pol(A(:,1),A(:,2))];
ang  = [ang ;-pi];

ang  = sort(unique(ang)); 


options = optimset('Display','off');
options = optimset(options,'LargeScale','off');

exitflag = 1;

 
   for i =1:length(ang)-1
      [c1,c2] = pol2cart((ang(i)+ang(i+1))/2,1);
      c = -[c1;c2];
      
      if i == 1      
         [X3(:,i),f,exitflag] = linprog(c,A,b,[],[],[],[],[],options);        
      else
         [X3(:,i),f,exitflag] = linprog(c,A,b,[],[],[],[],X3(:,i-1),options);
      end

      if and(i==1,exitflag ~=1)
         break
      end    
   end
   area  = polyarea(X3(1,:),X3(2,:));

