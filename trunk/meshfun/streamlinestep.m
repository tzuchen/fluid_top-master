function [x,w] = streamlinestep(x,u,ds,cdinfo)
% Given the velocity field as a vector, a location x and 
% a step size ds, this function finds the next x along 
% the streamline.
%
%  Tzu-Chen Liang 11/15/2005

  dim     = cdinfo.dim;
  nofgrid = cell2mat(cdinfo.nofgrid);
  w       = cell2mat(cdinfo.coord);
  gsize   = cell2mat({w.gsize})';
  binstr  = dec2bin(0:7);
  empstr  = ['        ']';
  bin     = str2num([binstr(:,1),empstr,binstr(:,2),empstr,binstr(:,3)])';  

 % If x has exceeded the boundry, we don't need to calculate it!
 xoffind=[];
 for i = 1:dim
       xoffind =unique([xoffind, union(find(x(i,:)<=cdinfo.coord{i}.l),find(x(i,:)>=cdinfo.coord{i}.r))]);
 end
 xind  = setdiff(1:size(x,2),xoffind);
 xr    = x(:,xind);
 xoff  = x(:,xoffind);


  lengthx = size(xr,2);

  % we need a larger w, because the ijk index could be imax+1,jmax+1,kmax+1
  gridl   = sum(nofgrid(1:dim))+nofgrid(dim)+1; 
  w       = sparse(lengthx,gridl); 

  for i = 1:dim
   
      [ijk,rem] = loc2ijk(xr,i,cdinfo);
      rem       = inv(diag(gsize))*rem;   
      vcijk     = kron(ijk,ones(1,8)) + kron(ones(1,lengthx),bin);
      vcind     = ijk2ind(vcijk,i,cdinfo);
      origin    = sum(nofgrid(1:(i-1)));
      wind      = origin+vcind;


      %[v,wt]    = lininterp3d(rem,reshape(u(wind),8,lengthx));

      [v,wt]    = lininterp3d(rem,reshape(soleval(u,wind,i,cdinfo),8,lengthx));
      xr(i,:)   = xr(i,:) + v*ds; 
      iind      = reshape(ones(8,1)*[1:lengthx],8*lengthx,1); 
      wi        = sparse(iind,wind,wt,lengthx,gridl);
      w         = w + wi;

  end

 w = w(:,1:sum(nofgrid(1:dim)));
 x(:,xind) = xr;
 x(:,xoffind) = xoff;

