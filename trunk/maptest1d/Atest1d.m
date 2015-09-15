function [x,mixnorm1,mixnorm2,mixnorm3,X]=Atest1d(A,ic,nit,X,showpic)

if nargin<2
   ic = 'sin';
end  
if nargin <3
   nit = 100;
end
if nargin <5
 showpic = 1;
end


n  = size(A,1);
dx = 1/n;
xlist = dx/2:dx:1-dx/2;

%%%%%%%%%%%%%%%5
% this is only true for logistic map
if nargin<4
 X = [];
end
if isempty(X)
 X = dx*(1/pi./(xlist.*(1-xlist)).^0.5)';
end
%%%%%%%%%%%%%%%


defaultic = (cos(2*pi*xlist)'+1)/2;
x0        = defaultic;
if isa(ic,'char')
   if ic == 'sinx'
       x0 = defaultic;
   end
   if ic == 'sing'
       x0 = 0*defaultic;
       x0(1:1) = 1;
   end  
   if ic == 'squx'
       x0 = 0*defaultic;
       x0(1:n/2) = 1;
   end   

   if ic == 'line'
       x0 = xlist';      
   end  

end
x0 = x0/sum(x0)*length(x0);

x=  x0;
mixnorm1 = [];
mixnorm2 = [];
mixnorm3 = [];


   y = logisticmap(xlist,'u')';

for i = 1:nit
 
    x = A*x; 
    p = x./y;
    meanp = sum(p)/n;
   % h = plot(xlist,x./y);
  if showpic == 1
    h = plot(xlist,x);
    axis([0 1 0 10]);
   
    M(i) = getframe;
  end
    xave= sum(x)/n;

    %1norm var (x as distribution)
    mixnorm1 = [mixnorm1 sum(abs(p-meanp)/n)];
    %TV        (x as distribution)
    mixnorm2 = [mixnorm2 sum(abs(x-y))];
    %func      (x as function)
    %mixnorm3 = [mixnorm3 1/2*sum(abs(x-sum(x)/n)/n)];
    %func      (x as function)
    mixnorm3 = [mixnorm3 1/2*sum(abs(x-sum(x)/n).*y/n)];   

end
%close all
mixnorm1 = mixnorm1./mixnorm1(1);
mixnorm2 = mixnorm2./mixnorm2(1);
%movie(M)
