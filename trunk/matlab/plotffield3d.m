function meshstruc  = plotffield(meshstruc,plotmat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coordinfo = meshstruc.coordinfo;
n         = meshstruc.n;
np        = meshstruc.np;
npt       = meshstruc.npt;
ngd       = meshstruc.ngd;
ugd       = meshstruc.ugd;
vgd       = meshstruc.vgd;
wgd       = meshstruc.wgd;
pgd       = meshstruc.pgd;
umat      = meshstruc.umat;
vmat      = meshstruc.vmat;
wmat      = meshstruc.wmat;


alphamat  = meshstruc.alphamat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uxcoord   = coordinfo.uxcoord;
uycoord   = coordinfo.uycoord;
uzcoord   = coordinfo.uzcoord;

vxcoord   = coordinfo.vxcoord;
vycoord   = coordinfo.vycoord;
vzcoord   = coordinfo.vzcoord;

wxcoord   = coordinfo.wxcoord;
wycoord   = coordinfo.wycoord;
wzcoord   = coordinfo.wzcoord;

pxcoord   = coordinfo.pxcoord;
pycoord   = coordinfo.pycoord;
pzcoord   = coordinfo.pzcoord;




figure
hold on
    
minx = min(uxcoord);
maxx = max(uxcoord);
miny = min(vycoord);
maxy = max(vycoord);
minz = min(wzcoord);
maxz = max(wzcoord);
    
     
     
%     for i = 1 : npt
%       h = line ( [uxcoord(i) uxcoord(i)], [ miny maxy ] );       
%       set(h, 'linewidth',0.5);
%     end
%     for j = 1 : npt
%       h = line ( [ minx maxx ], [vycoord(j) vycoord(j)] );  
%       set(h, 'linewidth',0.5,'color','g');
%     end
     
%      
%     for npcount = 1:np
%         h = plot3( [ pgd{npcount}.x pgd{npcount}.x+3*pgd{npcount}.u ],...
%                    [ pgd{npcount}.y pgd{npcount}.y+3*pgd{npcount}.v ],...
%                    [ pgd{npcount}.z pgd{npcount}.z+3*pgd{npcount}.w ],...
%                      'g');
%         set(h, 'linewidth',1);
%     end    
    
    
    pumat = zeros(ngd,ngd,ngd);
    pvmat = zeros(ngd,ngd,ngd);
    pwmat = zeros(ngd,ngd,ngd);
    for npcount = 1:np
        i = pgd{npcount}.i;
        j = pgd{npcount}.j;
        k = pgd{npcount}.k;
        pumat(i,j,k) = pgd{npcount}.u ;
        pvmat(i,j,k) = pgd{npcount}.v ;
        pwmat(i,j,k) = pgd{npcount}.w ;
        
        
    end      
    
    for ncount = 1:n
        if ugd{ncount}.c~=0
        h = plot3( [ ugd{ncount}.x ugd{ncount}.x+ugd{ncount}.c/10 ],...
                  [ ugd{ncount}.y ugd{ncount}.y                  ],...
                  [ ugd{ncount}.z ugd{ncount}.z                  ],...
                  'r');
        set(h, 'linewidth',1);
        end
        if vgd{ncount}.c~=0
        h = plot3( [ vgd{ncount}.x vgd{ncount}.x                  ],...
                  [ vgd{ncount}.y vgd{ncount}.y+vgd{ncount}.c/10 ],...
                  [ vgd{ncount}.z vgd{ncount}.z                  ],...
                  'r');
        set(h, 'linewidth',1);
        end
        if wgd{ncount}.c~=0
        h = plot3( [ wgd{ncount}.x wgd{ncount}.x                  ],...
                  [ wgd{ncount}.y wgd{ncount}.y                  ],...
                  [ wgd{ncount}.z wgd{ncount}.z+wgd{ncount}.c/10 ],...
                  'r');
        set(h, 'linewidth',1);
        end
    
    end    
    
    
    
%     if plotmat ==1
%         for ncount = 1:n 
%              if ugd{ncount}.alphax> 0.5
%                 h=plot(ugd{ncount}.x,ugd{ncount}.y,'or');
%                 set(h,'linewidth',2);
%              end       
%         end
%         for ncount = 1:n
%              if vgd{ncount}.alphay> 0.5
%                 h=plot(vgd{ncount}.x,vgd{ncount}.y,'or');
%                 set(h,'linewidth',2);
%              end       
%         end
%     end 
%     
%     
%     if plotmat == 2
%        figure
%        surf(alphamat')
% 
%     end
    
    
    pumat = permute(pumat,[2,1,3]);
    pvmat = permute(pvmat,[2,1,3]);
    pwmat = permute(pwmat,[2,1,3]);
    
     size(pumat)
     [X,Y,Z] = meshgrid(pxcoord,pycoord,pzcoord);
     [Yp,Zp] = meshgrid(pycoord,pzcoord);
     Xp = 0*Yp+0.1;
    % streamline(X,Y,Z,pumat,pvmat,pwmat,Xp,Yp,Zp)
     streamline(X,Y,Z,pumat,pvmat,pwmat,1*0.4*ones(9,1), 0.5*ones(9,1),0.1:0.1:0.9)
      streamline(X,Y,Z,-pumat,-pvmat,-pwmat,1*0.4*ones(9,1), 0.5*ones(9,1),0.1:0.1:0.9)
      
      streamline(X,Y,Z,pumat,pvmat,pwmat,1*0.4*ones(9,1),0.1:0.1:0.9, 0.5*ones(9,1))
      streamline(X,Y,Z,-pumat,-pvmat,-pwmat,1*0.4*ones(9,1),0.1:0.1:0.9, 0.5*ones(9,1))
      
       streamline(X,Y,Z,pumat,pvmat,pwmat,1*0.4*ones(9,1),0.1:0.1:0.9, 0.4*ones(9,1))
       streamline(X,Y,Z,-pumat,-pvmat,-pwmat,1*0.4*ones(9,1),0.1:0.1:0.9, 0.4*ones(9,1))
             streamline(X,Y,Z,pumat,pvmat,pwmat,1*0.4*ones(9,1),0.1:0.1:0.9, 0.6*ones(9,1))
       streamline(X,Y,Z,-pumat,-pvmat,-pwmat,1*0.4*ones(9,1),0.1:0.1:0.9, 0.6*ones(9,1))
%     contour(X,Y,alphamat');

     %p = patch(isosurface(X,Y,Z,permute(alphamat,[2,1,3]),[10]))
      p = patch(isosurface(X,Y,Z,permute(alphamat,[2,1,3]),[10]))
     isonormals(X,Y,Z,permute(alphamat,[2,1,3]), p)
     set(p, 'FaceColor', 'red', 'EdgeColor', 'none');
%      
%        
%     
       % daspect([1 1 1]); axis tight; 
        colormap(prism(28))
       % camup([1 0 0 ]); campos([25 -55 5]) 
        camlight; lighting phong
     
     
    axis equal
    
    axis([0 1 0 1 0 1])
    grid on
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

meshstruc.pumat = pumat;
meshstruc.pvmat = pvmat;
meshstruc.pwmat = pwmat;

    