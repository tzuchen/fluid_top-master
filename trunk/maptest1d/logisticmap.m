function y = logisticmap(x,type)
% Logistic Map
%  Type has following options:
%  empty : do tent map
%  'd' : y = S'(x),  the derivative of tent map at x 
%  'i' : y = S^-1(x), y has two rows..
%  'u' : the invariant distribution of S on x 
%
% Tzu-Chen Liang 5-18-2007 

c= 4;
% make x a row vector
if size(x,2)==1
  x = x';
end
 
if nargin==1
    y = c.*x.*(1-x);
else
   if type == 'i'  %  S^-1
     y = [(1-sqrt(1-4/c.*x))/2; (1+sqrt(1-4/c.*x))/2];              
   end
   if type == 'd'  % the derivative of S
     y = c*(1-2*x);  
   end
   if type == 'u'
       if c~=4
         sprintf('warning: this optin only valid for c=4')
       end      
       y = 1./pi./sqrt(x.*(1-x));

   end   

end
