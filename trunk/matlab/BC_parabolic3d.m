% this function add B.C to meshdata
% this is a test BC function for periodic solver
function meshstruc  =  BC_parabolic(meshstruc,info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ngd        = meshstruc.ngd;
 npt        = meshstruc.npt;
 dx         = meshstruc.dx;
 dy         = meshstruc.dy;
 dz         = meshstruc.dz;
 n          = meshstruc.n;
 np         = meshstruc.np;
 f          = meshstruc.f;
 ugd        = meshstruc.ugd;;
 vgd        = meshstruc.vgd;
 wgd        = meshstruc.wgd;
 pgd        = meshstruc.pgd;
 uindexmat  = meshstruc.uindexmat;
 vindexmat  = meshstruc.vindexmat;
 windexmat  = meshstruc.windexmat;
 pindexmat  = meshstruc.pindexmat;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B.C here (modify here to modify BC)

 for ncount = 1:n
        if ugd{ncount}.x ==1                         
              ugd{ncount}.bc   = 1;              
        end
        if or(vgd{ncount}.y ==0,vgd{ncount}.y ==1)                         
              vgd{ncount}.bc   = 1;            
        end
        if or(wgd{ncount}.z ==0,wgd{ncount}.z ==1)                         
              wgd{ncount}.bc   = 1;            
        end                

 end


 
 
 

 objtype  = 'linear';


 % objective function setting
 % Now it is set in ObjectiveSetting.m
  cx  = zeros(n,1);
  cy  = zeros(n,1);
  cz  = zeros(n,1);
     
  c = [cx;cy;cz];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% design parameter (alpha) setting  
alphalist = [];
for npcount = 1:np
    if or(pgd{npcount}.z < 0.2,pgd{npcount}.z > 0.8)
         alphalist = [ alphalist npcount ];
        
    end
end
alphalist = sort(alphalist);
nalpha = length(alphalist); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bodyforcex  =   1;
bodyforcey  =   0;
bodyforcez  =   0;

f         =   [ bodyforcex * ones(n,1); 
                bodyforcey * ones(n,1);
                bodyforcez * ones(n,1) ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalize the BC
bcnum = 0;
bcsum = 0;
 for ncount = 1:n
   if and(ugd{ncount}.bc==1,ugd{ncount}.u~=0) 
       bcnum = bcnum + 1;
       if ugd{ncount}.x ==1 
            bcsum = bcsum + ugd{ncount}.u;      
       else
            bcsum = bcsum - ugd{ncount}.u;
       end    
   end
 end
 
  for ncount = 1:n
   if and(vgd{ncount}.bc==1,vgd{ncount}.v~=0) 
       bcnum = bcnum + 1;
       if vgd{ncount}.y ==1 
            bcsum = bcsum + vgd{ncount}.v;      
       else
            bcsum = bcsum - vgd{ncount}.v;
       end    
   end
  end
  
  for ncount = 1:n
   if and(wgd{ncount}.bc==1,wgd{ncount}.w~=0) 
       bcnum = bcnum + 1;
       if wgd{ncount}.z ==1 
            bcsum = bcsum + wgd{ncount}.w;      
       else
            bcsum = bcsum - wgd{ncount}.w;
       end    
   end
  end
  
  
  
  
 if bcnum>0 
    bcaverage = bcsum/bcnum;
 else
    bcaverage = 0;
 end

 
 for ncount = 1:n
    if and(ugd{ncount}.bc==1,ugd{ncount}.u~=0)
        ugd{ncount}.u = ugd{ncount}.u - bcaverage ;
    end
 end
 for ncount = 1:n
    if and(vgd{ncount}.bc==1,vgd{ncount}.v~=0)
        vgd{ncount}.v = vgd{ncount}.v - bcaverage ;
    end
 end

 for ncount = 1:n
    if and(wgd{ncount}.bc==1,wgd{ncount}.w~=0)
        wgd{ncount}.w = wgd{ncount}.w - bcaverage ;
    end
 end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ncount = 1:n
    ugd{ncount}.c = c(ncount);
    vgd{ncount}.c = c(ncount+n);  
    wgd{ncount}.c = c(ncount+2*n);  
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 meshstruc.ngd       = ngd;
 meshstruc.npt       = npt;
 meshstruc.n         = n;
 meshstruc.nalpha    = nalpha;
 meshstruc.np        = np;
 meshstruc.f         = f; 
 meshstruc.ugd       = ugd;
 meshstruc.vgd       = vgd;
 meshstruc.wgd       = wgd;
 meshstruc.pgd       = pgd;
 meshstruc.c         = c;
 meshstruc.objtype   = objtype;
 meshstruc.alphalist = alphalist;






