function y = halfsinmap(x,type)
% Tent Map
%  Type has following options:
%  empty : do tent map
%  'd' : y = S'(x),  the derivative of tent map at x 
%  'i' : y = S^-1(x), y has two rows..
%  'u' : the invariant distribution of S on x 
%
% Tzu-Chen Liang 5-21-2007 

% make x a row vector
if size(x,2)==1
  x = x';
end
 
if nargin==1
    y = sin(pi*x);
else
   if type == 'i'  %  S^-1
     y = [ asin(x)/pi; 1-asin(x)/pi];              
   end
   if type == 'd'  % the derivative of S
     y = pi*cos(pi*x);  
   end
   if type == 'u'
     y = wfunS(x,8,1,@halfsinmap);       
   end
end
