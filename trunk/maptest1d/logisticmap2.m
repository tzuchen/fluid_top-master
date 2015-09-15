function y = logisticmap2(x,type)
% Logistic Map
%  Type has following options:
%  empty : do tent map
%  'd' : y = S'(x),  the derivative of tent map at x 
%  'i' : y = S^-1(x), y has two rows..
%  'u' : the invariant distribution of S on x 
%
% Tzu-Chen Liang 5-18-2007 

c= 16;
% make x a row vector
if size(x,2)==1
  x = x';
end
 
if nargin==1
    y = c.*x.^2.*(1-x).^2;
else
   if type == 'i'  %  S^-1
     y = [(1-sqrt(1-sqrt(x)))/2; (1+sqrt(1-sqrt(x)))/2];              
   end
   if type == 'd'  % the derivative of S
     y = 2*c*x.*(1-3*x+2*x.^2);  
   end
   if type == 'u'   
       y = wfunS(x,8,1,@logisticmap2);
   end   

end
