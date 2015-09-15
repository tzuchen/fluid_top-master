function meshstruc = generateoperator3ds(meshstruc,info)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ngd        = meshstruc.ngd;
 npt        = meshstruc.npt;
 dx         = meshstruc.dx;
 dy         = meshstruc.dy;
 dz         = meshstruc.dz;
 n          = meshstruc.n;
 np         = meshstruc.np;
 f          = meshstruc.f;
 c          = meshstruc.c;
 ugd        = meshstruc.ugd;
 vgd        = meshstruc.vgd;
 wgd        = meshstruc.wgd;
 pgd        = meshstruc.pgd;
 uindexmat  = meshstruc.uindexmat;
 vindexmat  = meshstruc.vindexmat;
 windexmat  = meshstruc.windexmat;
 pindexmat  = meshstruc.pindexmat;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vbclist       = [];
ubclist       = [];
wbclist       = [];
vbcvalue      = [];
ubcvalue      = [];
wbcvalue      = [];

for ncount = 1:npt*ngd*ngd
        if  ugd{ncount}.bc ==1
            ubclist   =   [ ubclist    ncount ];
            ubcvalue =    [ ubcvalue   ugd{ncount}.u ];
        end
end
for ncount = 1:ngd*npt*ngd
        if  vgd{ncount}.bc ==1
            vbclist   =   [ vbclist    ncount ];
            vbcvalue =    [ vbcvalue   vgd{ncount}.v ];
        end
end
for ncount = 1:ngd*ngd*npt
        if  wgd{ncount}.bc ==1
            wbclist   =   [ wbclist    ncount ];
            wbcvalue =    [ wbcvalue   wgd{ncount}.w ];
        end
end



[ubclist , ind ]  =  sort(ubclist);
ubcvalue          =  ubcvalue(ind);
usollist          =  setdiff(1:npt*ngd*ngd  , ubclist);
nur               =  length(usollist);

[vbclist , ind ]  =  sort(vbclist);
vbcvalue          =  vbcvalue(ind);
vsollist          =  setdiff(1:ngd*npt*ngd  , vbclist);
nvr               =  length(vsollist);

[wbclist , ind ]  =  sort(wbclist);
wbcvalue          =  wbcvalue(ind);
wsollist          =  setdiff(1:ngd*ngd*npt  , wbclist);
nwr               =  length(wsollist);


bclist            = [ ubclist   (n+vbclist) (n+n+wbclist)];
bcvalue           = [ ubcvalue  vbcvalue wbcvalue];

sollist           = setdiff(1:3*n, bclist);
nr                = nur+nvr+nwr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% form Lr and lr
L      = sparse(3*n, 3*n);

    for ncount = 1 : n
        i = ugd{ncount}.i;
        j = ugd{ncount}.j;
        k = ugd{ncount}.k;
        
        L(ncount,ncount) = -2/dx^2 - 2/dy^2 - 2/dz^2 ;
        if i > 1
            L(ncount,uindexmat(i-1,j  ,k  )) = 1/dx^2;
        end
        if i < npt
            L(ncount,uindexmat(i+1,j  ,k  )) = 1/dx^2;
        end
        if j > 1
            L(ncount,uindexmat(i  ,j-1,k  )) = 1/dy^2;
        end
        if j < ngd
            L(ncount,uindexmat(i  ,j+1,k  )) = 1/dy^2; 
        end     
        
        if k > 1
            L(ncount,uindexmat(i  ,j  ,k-1)) = 1/dz^2;
        end
        if k < ngd
            L(ncount,uindexmat(i  ,j  ,k+1)) = 1/dz^2; 
        end    
            
        
        
        if j == 1  
             %L(ncount,ncount) = -2/dx^2 - 4/dy^2;
             L(ncount,uindexmat(i ,j+1,k)) = 1/dy^2;%4/3/dy^2; 
        end
        if j == ngd
             %L(ncount,ncount) = -2/dx^2 - 4/dy^2;
             L(ncount,uindexmat(i ,j-1,k)) =1/dy^2;% 4/3/dy^2; 
        end
        
        if k == 1  
             %L(ncount,ncount) = -2/dx^2 - 4/dz^2;
             L(ncount,uindexmat(i ,j ,k+1)) = 1/dz^2;%4/3/dz^2; 
        end
        if k == ngd
             %L(ncount,ncount) = -2/dx^2 - 4/dz^2;
             L(ncount,uindexmat(i ,j ,k-1)) = 1/dz^2;%4/3/dz^2; 
        end
        
        
        
        
        % period condition
        if info.period == 1
            if i == 1
                L(ncount,uindexmat(npt-1, j, k)) = 1/dx^2;
            end
            if i == npt-1
                L(ncount,uindexmat(1, j, k))     = 1/dx^2;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for ncount = n+1 : 2*n
        i = vgd{ncount-n}.i;
        j = vgd{ncount-n}.j;
        k = vgd{ncount-n}.k;
        
        L(ncount,ncount) = -2/dx^2 - 2/dy^2 - 2/dz^2;
        if i > 1
            L(ncount,vindexmat(i-1,j  ,k  )) = 1/dx^2;
        end
        if i < ngd
            L(ncount,vindexmat(i+1,j  ,k  )) = 1/dx^2;
        end
        if j > 1
            L(ncount,vindexmat(i  ,j-1,k  )) = 1/dy^2;
        end
        if j < npt
            L(ncount,vindexmat(i  ,j+1,k  )) = 1/dy^2; 
        end   
        if k > 1
            L(ncount,vindexmat(i  ,j  ,k-1)) = 1/dz^2;
        end
        if k < ngd
            L(ncount,vindexmat(i  ,j  ,k+1)) = 1/dz^2; 
        end   
        
        
        
        
        if info.period ~= 1;
            if i == 1  
                 %L(ncount,ncount) = -2/dy^2 - 4/dx^2;
                 L(ncount,vindexmat(i+1 ,j  ,k  )) = 1/dz^2;%4/3/dy^2; 
            end
            if i == ngd
                 %L(ncount,ncount) = -2/dy^2 - 4/dx^2;
                 L(ncount,vindexmat(i-1 ,j  ,k  )) = 1/dz^2;%4/3/dy^2; 
            end

            if k == 1  
                 %L(ncount,ncount) = -2/dy^2 - 4/dz^2;
                 L(ncount,vindexmat(i  ,j  ,k+1)) = 1/dz^2;%4/3/dz^2; 
            end
            if k == ngd
                 %L(ncount,ncount) = -2/dy^2 - 4/dz^2;
                 L(ncount,vindexmat(i  ,j  ,k-1)) = 1/dz^2;%4/3/dz^2; 
            end
        end
        
        
        
        if info.period == 1  
        % period condition  
            if i == 1 
                L(ncount,vindexmat(ngd, j, k)) = 1/dx^2;
            end
            if i == ngd
                L(ncount,vindexmat(1, j, k)) = 1/dx^2;
            end
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for ncount = 2*n+1 : 3*n
        i = wgd{ncount-2*n}.i;
        j = wgd{ncount-2*n}.j;
        k = wgd{ncount-2*n}.k;
        
        L(ncount,ncount) = -2/dx^2 - 2/dy^2 - 2/dz^2;
        if i > 1
            L(ncount,windexmat(i-1,j  ,k  )) = 1/dx^2;
        end
        if i < ngd
            L(ncount,windexmat(i+1,j  ,k  )) = 1/dx^2;
        end
        if j > 1
            L(ncount,windexmat(i  ,j-1,k  )) = 1/dy^2;
        end
        if j < ngd
            L(ncount,windexmat(i  ,j+1,k  )) = 1/dy^2; 
        end   
        if k > 1
            L(ncount,windexmat(i  ,j  ,k-1)) = 1/dz^2;
        end
        if k < npt
            L(ncount,windexmat(i  ,j  ,k+1)) = 1/dz^2; 
        end   
        
        
        
        
        if info.period ~= 1;
            if i == 1  
                 %L(ncount,ncount) = -2/dz^2 - 4/dx^2;
                 L(ncount,windexmat(i+1 ,j  ,k  )) = 1/dy^2;%4/3/dy^2; 
            end
            if i == ngd
                 %L(ncount,ncount) = -2/dz^2 - 4/dx^2;
                 L(ncount,windexmat(i-1 ,j  ,k  )) = 1/dy^2;%4/3/dy^2; 
            end

            if j == 1  
                 %L(ncount,ncount) = -2/dz^2 - 4/dy^2;
                 L(ncount,windexmat(i  ,j+1,k  )) = 1/dz^2;%4/3/dz^2; 
            end
            if j == ngd
                 %L(ncount,ncount) = -2/dz^2 - 4/dy^2;
                 L(ncount,windexmat(i  ,j-1,k  )) = 1/dz^2;%4/3/dz^2; 
            end
        end
        
        
        
        if info.period == 1  
        % period condition  
            if i == 1 
                L(ncount,windexmat(ngd, j , k)) = 1/dx^2;
            end
            if i == ngd
                L(ncount,windexmat(1, j ,k)) = 1/dx^2;
            end
        end
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
Lr                             =  L(sollist,sollist);
lr                             =  L(sollist, bclist ) * bcvalue';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% form Dr and dr, Gr 
D = sparse(np, 3*n);
    for npcount = 1 :ngd*ngd*ngd
        
        i           = pgd{npcount}.i;
        j           = pgd{npcount}.j;
        k           = pgd{npcount}.k;
        
        D(pindexmat(i,j,k),uindexmat(i+1,j  ,k  )) =   1/dx;
        D(pindexmat(i,j,k),uindexmat(i  ,j  ,k  )) =  -1/dx;        
        D(pindexmat(i,j,k),vindexmat(i  ,j+1,k  )) =   1/dy;
        D(pindexmat(i,j,k),vindexmat(i  ,j  ,k  )) =  -1/dy;
        D(pindexmat(i,j,k),windexmat(i  ,j  ,k+1)) =   1/dz;
        D(pindexmat(i,j,k),windexmat(i  ,j  ,k  )) =  -1/dz;
        
        if info.period == 1 
        % period condition
            if i == ngd
                  D(pindexmat(i,j,k),uindexmat(1,j,k  ))   =   1/dx;
                  D(pindexmat(i,j,k),uindexmat(i+1,j,k  )) =   0;
            end
        end
        
    end
   D =  -D;

   Dr        =  D(:, sollist);
   dr        =  D(:, bclist ) * bcvalue';
   G         =  D';
   Gr        =  Dr';
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fr        = f(sollist);
  cr        =  [ c(sollist); zeros(np,1) ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P is the matrix maps P space to UV space
P = sparse(3*n,np);
for ncount = 1:n
     i = ugd{ncount}.i;
     j = ugd{ncount}.j;
     k = ugd{ncount}.k;
    
     if i ==1
          % periodic condition is considered
          P(ncount,pindexmat(i  ,j,k)) = 1/2;
          P(ncount,pindexmat(ngd,j,k)) = 1/2;           
     elseif i == npt
          % Do nothing         
     else 
          P(ncount,pindexmat(i-1,j,k)) = 1/2;
          P(ncount,pindexmat(i  ,j,k)) = 1/2;
     end     
end

for ncount = 1:n
     i = vgd{ncount}.i;
     j = vgd{ncount}.j;
     k = vgd{ncount}.k;
    
     if j ==1         
          P(ncount+n,pindexmat(i  ,j,k)) = 1;                     
     elseif j == npt
          P(ncount+n,pindexmat(i,ngd,k)) = 1;         
     else 
          P(ncount+n,pindexmat(i,j-1,k)) = 1/2;
          P(ncount+n,pindexmat(i,j  ,k)) = 1/2;
     end     
end


for ncount = 1:n
     i = wgd{ncount}.i;
     j = wgd{ncount}.j;
     k = wgd{ncount}.k;
    
     if k ==1         
          P(ncount+2*n,pindexmat(i  ,j,k)) = 1;                     
     elseif k == npt
          P(ncount+2*n,pindexmat(i,j,ngd)) = 1;         
     else 
          P(ncount+2*n,pindexmat(i,j  ,k-1)) = 1/2;
          P(ncount+2*n,pindexmat(i,j  ,k  )) = 1/2;
     end     
end

Pr = P(sollist,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alphavec = sparse(np,1);
alphak   = sparse(nr,1);
objlist  = [];

sol    = sparse(nr+np,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 meshstruc.fr         = fr; 
 % new data
 meshstruc.Lr        = Lr;
 meshstruc.lr        = lr;
 meshstruc.Dr        = Dr;
 meshstruc.dr        = dr;
 meshstruc.Gr        = Gr;
 meshstruc.Pr        = Pr;
 meshstruc.nr        = nr; 
 meshstruc.nur       = nur;
 meshstruc.nvr       = nvr;
 meshstruc.nwr       = nwr;
 meshstruc.sollist   = sollist;
 meshstruc.bclist    = bclist;
 meshstruc.usollist  = usollist;
 meshstruc.vsollist  = vsollist;
 meshstruc.wsollist  = wsollist;
 meshstruc.bcvalue   = bcvalue;
 meshstruc.ubcvalue  = ubcvalue;
 meshstruc.vbcvalue  = vbcvalue;
 meshstruc.wbcvalue  = wbcvalue;
 meshstruc.alphak    = alphak;
 meshstruc.alphavec  = alphavec;
 meshstruc.sol       = sol;
 meshstruc.objlist   = [];
 meshstruc.trueobjlist = [];
 meshstruc.cr        = cr;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






  
