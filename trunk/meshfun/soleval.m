function x = soleval(u,ind,dir,cdinfo)
% This function returns x = u(ind)
% However, it checks if each ind is inside the correct range
% of u, if not, the corresponding x is zero


nofgrid = cell2mat(cdinfo.nofgrid);
Istart = sum(nofgrid(1:(dir-1)))+1;
Iend    = sum(nofgrid(1:dir));


  zeroset = union(find(ind<Istart),find(ind>Iend));
  x = u(ind);
  x(zeroset) = 0;

