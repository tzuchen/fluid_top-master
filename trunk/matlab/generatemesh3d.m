% generate mesh
function   meshstruc = generatemesh(ngd,info)

npt    =    ngd+1;

dx     =    1/ngd;
dy     =    1/ngd;
dz     =    1/ngd;

mu     =    1;
bodyforcex  =  0;
bodyforcey  =  0;
bodyforcez  =  0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n         =   npt*ngd*ngd;
np        =   ngd*ngd*ngd;

f         =   [ bodyforcex * ones(n,1); 
                bodyforcey * ones(n,1);
                bodyforcez * ones(n,1) ];

uxcoord    =  0    : dx : 1;
uycoord    =  dy/2 : dy : 1-dy/2;
uzcoord    =  dz/2 : dz : 1-dz/2;

vxcoord    =  dx/2 : dx : 1-dx/2;
vycoord    =  0    : dy : 1;
vzcoord    =  dz/2 : dz : 1-dz/2;
 
wxcoord    =  dx/2 : dx : 1-dx/2;
wycoord    =  dy/2 : dy : 1-dy/2;
wzcoord    =  0    : dz : 1;


pxcoord   =  dx/2  :dx  : 1-dx/2;
pycoord   =  dy/2  :dy  : 1-dy/2;
pzcoord   =  dz/2  :dz  : 1-dz/2;


coordinfo.uxcoord = uxcoord;
coordinfo.uycoord = uycoord;
coordinfo.uzcoord = uzcoord;

coordinfo.vxcoord = vxcoord;
coordinfo.vycoord = vycoord;
coordinfo.vzcoord = vzcoord;

coordinfo.wxcoord = wxcoord;
coordinfo.wycoord = wycoord;
coordinfo.wzcoord = wzcoord;



coordinfo.pxcoord = pxcoord;
coordinfo.pycoord = pycoord;
coordinfo.pzcoord = pzcoord;




uindexmat  = zeros(npt,ngd,ngd);
vindexmat  = zeros(ngd,npt,ngd);
windexmat  = zeros(ngd,ngd,npt);

pindexmat  = zeros(ngd,ngd,ngd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncount   = 0;
for k = 1:ngd
    for j = 1:ngd
        for i = 1:npt
        
            ncount               = ncount + 1;
            uindexmat(i,j,k)     = ncount;
            ugd{ncount}.i        = i;
            ugd{ncount}.j        = j;
            ugd{ncount}.k        = k;
            ugd{ncount}.x      = uxcoord(i);
            ugd{ncount}.y      = uycoord(j);
            ugd{ncount}.z      = uzcoord(k);  
            ugd{ncount}.u      = 0;       
            ugd{ncount}.bc     = 0;
            ugd{ncount}.alphax = 0;
            ugd{ncount}.c      = 0;
        end
    end
end


imat = 


ncount   = 0;
for k = 1:ngd
    for j = 1:npt
        for i = 1:ngd
            ncount               = ncount + 1;
            vindexmat(i,j,k)     = n + ncount;
            vgd{ncount}.i        = i;
            vgd{ncount}.j        = j;
            vgd{ncount}.k        = k;
            vgd{ncount}.x        = vxcoord(i);
            vgd{ncount}.y        = vycoord(j);
            vgd{ncount}.z        = vzcoord(k);
            vgd{ncount}.v        = 0;       
            vgd{ncount}.bc       = 0;    
            vgd{ncount}.alphay   = 0;
            vgd{ncount}.c        = 0;
        end 
    end
end


ncount   = 0;
for k = 1:npt
    for j = 1:ngd
        for i = 1:ngd
            ncount               = ncount + 1;
            windexmat(i,j,k)     = n + n + ncount;
            wgd{ncount}.i        = i;
            wgd{ncount}.j        = j;
            wgd{ncount}.k        = k;
            wgd{ncount}.x        = wxcoord(i);
            wgd{ncount}.y        = wycoord(j);
            wgd{ncount}.z        = wzcoord(k);
            wgd{ncount}.w        = 0;       
            wgd{ncount}.bc       = 0;    
            wgd{ncount}.alphay   = 0;
            wgd{ncount}.c        = 0;
        end 
    end
end




npcount = 0;
for k = 1:ngd
    for j = 1:ngd
        for i = 1:ngd
            npcount           = npcount + 1;
            pindexmat(i,j,k)    = npcount;
            pgd{npcount}.i    = i;
            pgd{npcount}.j    = j;
            pgd{npcount}.k    = k;
            pgd{npcount}.x    = pxcoord(i);
            pgd{npcount}.y    = pycoord(j);
            pgd{npcount}.z    = pzcoord(k);
            pgd{npcount}.p    = 0;
            pgd{npcount}.c    = 0;   
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set velocity at boundary =0
for ncount = 1:npt*ngd
    if info.period~=1
        if or(ugd{ncount}.x ==0,ugd{ncount}.x ==1)                         
              ugd{ncount}.bc   = 1;              
        end
        if or(vgd{ncount}.y ==0,vgd{ncount}.y ==1)                         
              vgd{ncount}.bc   = 1;              
        end
        if or(wgd{ncount}.z ==0,wgd{ncount}.z ==1)                         
              wgd{ncount}.bc   = 1;              
        end
    end            
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alphamat = [];
adjvec   = [];
sol      = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 meshstruc.ngd       = ngd;
 meshstruc.npt       = npt;
 meshstruc.dx        = dx;
 meshstruc.dy        = dy;
 meshstruc.dz        = dz;
 meshstruc.n         = n;
 meshstruc.np        = np;
 meshstruc.f         = f; 
 meshstruc.mu        = mu;
 meshstruc.ugd       = ugd;
 meshstruc.vgd       = vgd;
 meshstruc.wgd       = wgd;
 meshstruc.pgd       = pgd;
 meshstruc.uindexmat = uindexmat;
 meshstruc.vindexmat = vindexmat;
 meshstruc.windexmat = windexmat;
 meshstruc.pindexmat = pindexmat;
 meshstruc.coordinfo = coordinfo;
 meshstruc.alphamat  = alphamat;
 meshstruc.sol       = sol;
 meshstruc.adjvec    = adjvec;
 






