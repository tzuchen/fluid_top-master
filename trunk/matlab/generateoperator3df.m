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

    [Lr,Dr,Gr,Pr,lr,dr]=mylaplacefun(ngd);
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fr        = f(sollist);
  cr        =  [ c(sollist); zeros(np,1) ];
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






  
