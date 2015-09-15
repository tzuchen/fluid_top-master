function y = tentmap(x,type)
% Tent Map
%  Type has following options:
%  empty : do tent map
%  'd' : y = S'(x),  the derivative of tent map at x 
%  'i' : y = S^-1(x), y has two rows..
%  'u' : the invariant distribution of S on x 
%
% Tzu-Chen Liang 5-18-2007 

% make x a row vector
%if size(x,2)==1
%  x = x';
%end
 
if nargin==1
    y = 1-2*abs(x-0.5);
else
   if type == 'i'  %  S^-1
     y = [ x/2; 1-x/2];              
   end
   if type == 'd'  % the derivative of S
     y = 4*((x<1/2)-1/2);  
   end
   if type == 'u'
     y = 0*x+1; 
   end
end

  

