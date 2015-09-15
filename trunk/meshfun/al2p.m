function al2pmat = al2p(alphaind,cdinfo)
% The function finds the map from alpha index
% to pressure grids
%
% Tzu-chen Liang 11/6/2005
  
  t0       = clock;
  nofpgrid = cdinfo.nofgrid{end};

  for setn = 1: size(alphaind,2)
     n =length(alphaind{setn});
     al2pmati{setn} = sparse(alphaind{setn},1:n,ones(n,1),nofpgrid,n);
  end

   al2pmat = cell2mat(al2pmati);

  if cdinfo.showtime==1
    disp(sprintf('Matlab: Time to generate al2p matrix : %f sec ',etime(clock,t0))); 
  end
