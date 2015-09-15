function be2almat = be2al(betadir,alphaind,cdinfo)
% This function returns the mapping from beta to alpha
%
% Tzu-Chen Liang 11/7/2005

t0   = clock;
n    = size(alphaind,2);
pdim = cdinfo.dim+1;


for setn = 1:n      
      indset      = ind2ijk(alphaind{setn},pdim,cdinfo);
      p           = length(indset);
      betaset     = unique(indset(betadir{setn},:));
      q           = length(betaset);        
      be2almati{setn} = kron(ones(q,1),speye(ceil(p/q)));
      perm        = [];

      [nouse, perm] = sort(indset(betadir{setn},:));
      r(perm)  = 1:p;

      be2almati{setn} = be2almati{setn}(r,:);

      r = [];
end


be2almat =  blkdiag(be2almati{:});

if cdinfo.showtime==1
  disp(sprintf('Matlab: Time to generate be2al matrix : %f sec ',etime(clock,t0))); 
end
