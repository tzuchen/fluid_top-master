function meshstruc  = beta2alpha2(meshstruc)

betavec = meshstruc.betavec;
np      = meshstruc.np;
pgd     = meshstruc.pgd;
Mind1   = meshstruc.pindexmat(:,:,1);
Mind2   = meshstruc.pindexmat(:,:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alphavec = zeros(np,1);


for pcount = 1:np
    z = pgd{pcount}.z;
  if z < 0.2  
    i = pgd{pcount}.i; 
    j = pgd{pcount}.j;
    s = betavec(Mind1(i,j)); 
    alphavec(pcount) = s;
  end
  if z > 0.8  
    i = pgd{pcount}.i; 
    j = pgd{pcount}.j;
    s = betavec(Mind2(i,j)); 
    alphavec(pcount) = s;
  end  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meshstruc.alphavec  =  alphavec;