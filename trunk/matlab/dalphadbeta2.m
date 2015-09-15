function meshstruc  = dalphadbeta(meshstruc)

pgd   = meshstruc.pgd;
np    = meshstruc.np;
ngd   = meshstruc.ngd;
Mind1 = meshstruc.pindexmat(:,:,1);
Mind2 = meshstruc.pindexmat(:,:,2);

if ~isfield(meshstruc,'betavec')
   meshstruc.betavec = zeros(ngd^2 *2+2,1); 
end
betavec = meshstruc.betavec;
nbeta   = length(betavec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M  = sparse(np,nbeta);

for pcount = 1:np
    if pgd{pcount}.z < 0.2;
       i = pgd{pcount}.i;
       j = pgd{pcount}.j;
       %k = pgd{pcount}.k;
       %z = pgd{pcount}.z;   
       %s = betavec(Mind1(i,j));
       M(pcount,Mind1(i,j)) =  1;
    end
    if pgd{pcount}.z > 0.8;
       i = pgd{pcount}.i;
       j = pgd{pcount}.j;
       %k = pgd{pcount}.k;
       %z = pgd{pcount}.z;
       %s = betavec(Mind2(i,j));
       M(pcount,Mind2(i,j)) = 1;
    end
    
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meshstruc.betavec = betavec;
meshstruc.nbeta   = nbeta;
meshstruc.M       = M;


