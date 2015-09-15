function [A,B,M,mdot] = map2markov(n,xpi,mapfunction,varargin)
% Given a map, this function returns the Markov matrix
% A and B /in R^{n*n} that approximate this map.
% 
% mapfunction has to be a function handle
% see standardmap.m for example 
% It has to be periodical and has period 1. 
%
% Tzu-Chen Liang 3-24-2006
%

 dx       = 1/n; 
 dy       = 1/n;
 x        = dx/2:dx:1-dx/2;
 y        = dy/2:dy:1-dy/2;

 [X,Y]    = meshgrid(x,y);
 xc       = reshape(X,n^2,1);
 yc       = reshape(Y,n^2,1);
 param    = varargin;
 [xce,yce]  = feval(mapfunction,xc,yc,0,1,param);

 xe       = 0:dx:1;
 ye       = 0:dy:1;
 [Xe,Ye]  = meshgrid(xe,ye);
 [Xe,Ye]  = feval(mapfunction,Xe,Ye,0,1,param);

 xe1      = reshape(Xe(1:end-1,1:end-1),n^2,1);
 ye1      = reshape(Ye(1:end-1,1:end-1),n^2,1);
 xe2      = reshape(Xe(2:end  ,1:end-1),n^2,1);
 ye2      = reshape(Ye(2:end  ,1:end-1),n^2,1);
 xe3      = reshape(Xe(1:end-1,2:end  ),n^2,1);
 ye3      = reshape(Ye(1:end-1,2:end  ),n^2,1);
 xe4      = reshape(Xe(2:end  ,2:end  ),n^2,1);
 ye4      = reshape(Ye(2:end  ,2:end  ),n^2,1);

 xl       =  min(min(ceil(xe1/dx),ceil(xe3/dx)),min(ceil(xe2/dx),ceil(xe4/dx)));
 xr       =  max(max(ceil(xe2/dx),ceil(xe4/dx)),max(ceil(xe1/dx),ceil(xe3/dx)));
 yd       =  min(min(ceil(ye3/dy),ceil(ye4/dy)),min(ceil(ye1/dy),ceil(ye2/dy)));
 yu       =  max(max(ceil(ye1/dy),ceil(ye2/dy)),max(ceil(ye3/dy),ceil(ye4/dy)));

nzind     = [];

for i = 1:n^2

       xlist     = mod(xl(i):xr(i),n);
       ylist     = mod(yd(i):yu(i),n);
       xlist(find(xlist==0))=n;
       ylist(find(ylist==0))=n;


       [Xi,Yi]   = meshgrid(xlist,ylist);
       xi        = reshape(Xi,1,prod(size(Xi)));
       yi        = reshape(Yi,1,prod(size(Yi)));
       nzind     = [ nzind ,((i-1)*n^2)+(xi-1)*n+yi];

end

 nzind   = unique(sort(nzind));
 nnzA    = length(nzind);
 [iind,jind] = ind2sub([n^2,n^2],nzind);


 R1      = sparse(iind,1:nnzA,ones(1,nnzA),n^2,nnzA);
 R2      = sparse(jind,1:nnzA,ones(1,nnzA),n^2,nnzA);
 R       = [R1 ; R2];
 b       = [xpi ; xpi];%1* ones(size(R,1),1);  
 f       = 0*ones(nnzA,1);  




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Mv = zeros(nnzA,1);
 for k = 1:nnzA
   i = iind(k);
   j = jind(k);
   Mv(k) = abs( dx-abs(xc(j)-xce(i)) )*abs( dy-abs(yc(j)-yce(i)) );
   
 end
   
  Mv=Mv/sum(Mv)*sum(xpi);
 

 Mp = sparse(n^2,n^2);
 Mp(nzind) = Mv;
  

 % Cx<d
  C  = [-speye(nnzA), -sparse(ones(nnzA,1)) ;
         speye(nnzA), -sparse(ones(nnzA,1))  ];
  d  =  [ -Mv; Mv ];
   
  R  = [R sparse(zeros(2*n^2,1))];
  f  = [f ;1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




 tic
    [xopt,fv,flag] = linprog(f,C,d,R,b,0*ones(nnzA,1),[]);
 toc
  flag
size(xopt)
size(R1)
  mdot           = R1*xopt(1:end-1);
  xopt           = xopt(1:end-1);

  Mvec           = sparse(n^4,1);
  Mvec(nzind)    = xopt;
  M              = reshape(Mvec,n^2,n^2);
  
  Avec           = sparse(n^4,1);  
  Bvec           = sparse(n^4,1); 

  for i = 1:nnzA
    k = nzind(i);
    Avec(k) = Mvec(k)/ mdot(iind(i));
    Bvec(k) = Mvec(k)/ mdot(jind(i));
  end
  A              = reshape(Avec,n^2,n^2);
  B              = reshape(Bvec,n^2,n^2)';
 
plot(Mv)
 hold on
  plot(M(nzind),'r')
 

  





