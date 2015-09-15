function [x,ind]  = ind2loc(i,dim,cdinfo)
% Given and a coordinate system cdinfo and an index i,
% return the location in space.
% Works for n-dimension
%    
%  Tzu-Chen Liang  11/6/2005


if isempty(i)
  x     = [];
  ind   = [];
else
  ind   = ind2ijk(i,dim,cdinfo);
  gsize = zeros(cdinfo.dim,1);

  for j = 1: cdinfo.dim
      gsize(j) = cdinfo.coord{j}.gsize;  
  end

  for j = 1: cdinfo.dim
      x(j,:) =gsize(j).*(ind(j,:)-0.5); 

  end

 if dim <= cdinfo.dim
    x(dim,:) = x(dim,:) + gsize(dim)/2;
    if cdinfo.coord{dim}.p==1
        x(dim,:) = x(dim,:) - gsize(dim);
    end
 end
  
end

