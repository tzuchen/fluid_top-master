function alphaind = alphaset(alphapos,cdinfo)
%  This function sets the alpha in the pressure grid points
%
%  Tzu-Chen Liang  11/6/2005
  
  t0   = clock;
  dim  = cdinfo.dim;
  pdim = dim+1; 
 

  nofalphaset = size(alphapos,2);
  pgrid  = ind2loc(1:cdinfo.nofgrid{pdim},pdim,cdinfo);

  for setn = 1: nofalphaset
      indexh = 1:cdinfo.nofgrid{pdim};
      for i = 1:dim         
          lb = alphapos{setn}(2*i-1);
          ub = alphapos{setn}(2*i  );
          indexl = indexh(find(lb<pgrid(i,indexh)));
          indexh = indexl(find(pgrid(i,indexl)<ub));
      end 
      alphaind{setn} = indexh;     
  end

if cdinfo.showtime==1
  disp(sprintf('Matlab: Time to set alpha : %f sec ',etime(clock,t0))); 
end
