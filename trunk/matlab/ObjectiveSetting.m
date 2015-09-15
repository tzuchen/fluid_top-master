function meshstruc  =  ObjectiveSetting(meshstruc,info)
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
 %uindexmat  = meshstruc.uindexmat;
 %vindexmat  = meshstruc.vindexmat;
 %windexmat  = meshstruc.windexmat;
 %pindexmat  = meshstruc.pindexmat;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  objtype  = 'linear';
 % objective function setting
  cx  = zeros(n,1);
  cy  = zeros(n,1);
  cz  = zeros(n,1);
%      for ncount = 1:n
%          if and(vgd{ncount}.x>0.4 ,vgd{ncount}.x<0.6)
%             if and(vgd{ncount}.y>0.4 ,vgd{ncount}.y<0.6)
%                if and(vgd{ncount}.z>0.5 ,wgd{ncount}.z<1) 
%                     cy(ncount) = 1;
%                end
%                 if and(vgd{ncount}.z>0.0 ,vgd{ncount}.z<0.5) 
%                      cy(ncount) = -1;
%                 end
% 
%             end
%         end     
%      end
      for ncount = 1:n
          if and(wgd{ncount}.x>0 ,wgd{ncount}.x<1)            
              if and(wgd{ncount}.z>0.4 ,wgd{ncount}.z<0.6) 
                if and(wgd{ncount}.y>0.7 ,wgd{ncount}.y<1)
                     cz(ncount) = 1;
                end
                if and(wgd{ncount}.y>0 ,wgd{ncount}.y<0.3) 
                      cz(ncount) = -1;
                end
 
             end
         end     
      end
     
%for ncount = 1:n
%      if  and(wgd{ncount}.z>0.4 ,wgd{ncount}.z<0.6) 
%          if and(wgd{ncount}.y>0 ,wgd{ncount}.y<0.1)
%                cz(ncount) = 1;
%          end
%          if and(wgd{ncount}.y>0.5 ,wgd{ncount}.y<0.7)
%                cz(ncount) = -1;
%          end
%          if and(wgd{ncount}.y>0.9 ,wgd{ncount}.y<1)
%                cz(ncount) = 1;
%          end
%          
%      end
%end

     
  c = [cx;cy;cz];
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ncount = 1:n
    ugd{ncount}.c = c(ncount);
    vgd{ncount}.c = c(ncount+n);  
    wgd{ncount}.c = c(ncount+2*n);     
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 meshstruc.ugd       = ugd;
 meshstruc.vgd       = vgd;
 meshstruc.wgd       = wgd;
 meshstruc.pgd       = pgd;
 meshstruc.c         = c;
 meshstruc.objtype   = objtype;

  
