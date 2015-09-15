 %
 % Tzu-Chen Liang   
 %
 % In this routine we generate the following operators for 3D
 % Stoke's flow using kronecker product. This is fast for huge 
 % problem. 
 % We assume nx = ny = nz = n here. 
 % 
 % Notice: All the operators need to be scaled by grid size 
 % before using! 
 % 
 % Lrf : Laplace Operator
 % Grf : Gradient Operator
 % Drf : Differential Operator
 % Prf : Pressure grid to velocity grid Operator
 %
 % lrf : The rhs terms of BC induced by Lrf
 % grf : The rhs terms of BC induced by Grf  
 %

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Information are given here.  
   n  = 4;
   dx = 1/n;
 % We assume it is periodical in x-direction
 % On the other 4 walls, BC of u,v and w need to be specified.  
 %
 %
 %
   n2  = n*n;
   n2r = n*(n-1);  

 % xy0u, xy1u
   xy0u = zeros(n2 ,1);
   xy1u = zeros(n2 ,1);
 
 % xy0v, xy1v
   xy0v = zeros(n2r,1);
   xy1v = zeros(n2r,1);
 
 % xy0w, xy1w
   xy0w = zeros(n2,1);
   xy1w = zeros(n2,1); 
 
 % zx0u, zx1u
   zx0u = zeros(n2,1);
   zx1u = zeros(n2,1);
 
 % zx0v, zx1v
   zx0v = zeros(n2,1);
   zx1v = zeros(n2,1);
 
 % zx0w, zx1w
   zx0w = zeros(n2r,1);
   zx1w = zeros(n2r,1);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % 1-D Operators 

   % L is a 1-D Laplace operator
   % Lp is the periodical version of it
   % Ls is a shrinked L
   % D is a 1-D Differential operator
   % Dp is the periodical version of it
   % Ds is a shrinked D
   % 
   % Notice: Lp and Dp are for x-direction, the periodical condition
   %         Ls and Ds are for the BCs. 

   I  	    = speye(n,n);
   E  	    = sparse(2:n,1:n-1,1,n,n);
   e  	    = ones(n,1);
   L        = spdiags([e -2*e e], -1:1, n, n);
    Lp       = L; 
   Lp(1,n)  = 1;
   Lp(n,1)  = 1;   
   Ls       = L(1:n-1,1:n-1);   

   D        = I-E';
   Dp       = D;
   Dp(n,1)  = -1;
   Ds       = D(:,2:n);
   
   Is       = speye(n-1,n-1);
   Ih       = I(:,1:n-1);
   Il       = I(:,2:n);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % All matrix operators

 % Laplace Operator
 %       x-dir               y-dir              z-dir
 L1   =  kron3(I , I , Lp) + kron3(I , L , I) + kron3(L , I , I);
 L2   =  kron3(I , Is, Lp) + kron3(I , Ls, I) + kron3(L , Is, I);
 L3   =  kron3(Is, I , Lp) + kron3(Is, L , I) + kron3(Ls, I,  I);
 Lrf  =  blkdiag(L1,L2,L3);
 
 
 % Differential Operator 
 D1   =  kron3(I , I  ,Dp);
 D2   =  kron3(I , Ds ,I );
 D3   =  kron3(Ds, I  ,I );
 Drf  =  [D1 D2 D3];
 
 % Gradient Operator
 Grf  = Drf';
 
 % P grid to UVW grid Operator
 Prf = abs(Drf')/2;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now deal with BCs 
  e0       = sparse(n,1);
  e0(1)    = 1;
  e1 	   = sparse(n,1);
  e1(n)    = 1; 
  e0s      = sparse(n-1,1);
  e0s(1)   = 1;
  e1s      = sparse(n-1,1);
  e1s(n-1) = 1;
  z0       = sparse(n,1);
  z0s      = sparse(n-1,1);
  z1       = sparse(n,1);
  z1s      = sparse(n-1,1);


  BCxy0u = kron3(e0 , I  , I );
  BCxy1u = kron3(e1 , I  , I );
  BCxy0v = kron3(e0 , Is , I ); 
  BCxy1v = kron3(e1 , Is , I ); 
  BCxy0w = kron3(e0s, I  , I ); 
  BCxy1w = kron3(e1s, I  , I );

  BCzx0u = kron3( I , e0 , I ); 
  BCzx1u = kron3( I , e1 , I ); 
  BCzx0v = kron3( I , e0s, I ); 
  BCzx1v = kron3( I , e1s, I ); 
  BCzx0w = kron3( Is, e0 , I ); 
  BCzx1w = kron3( Is, e1 , I );

  %                z   y    x
  lrfxy0u = BCxy0u*xy0u; 
  lrfxy1u = BCxy1u*xy1u; 
  lrfxy0v = BCxy0v*xy0v; 
  lrfxy1v = BCxy1v*xy1v; 
  lrfxy0w = BCxy0w*xy0w; 
  lrfxy1w = BCxy1w*xy1w;
 
  lrfzx0u = BCzx0u*zx0u; 
  lrfzx1u = BCzx1u*zx1u; 
  lrfzx0v = BCzx0v*zx0v; 
  lrfzx1v = BCzx1v*zx1v; 
  lrfzx0w = BCzx0w*zx0w; 
  lrfzx1w = BCzx1w*zx1w;
  
  lrf = [lrfxy0u + lrfxy1u + lrfzx0u + lrfzx1u;
         lrfxy0v + lrfxy1v + lrfzx0v + lrfzx1v;
         lrfxy0w + lrfxy1w + lrfzx0w + lrfzx1w ];
         
  % has problem here

  BCxy0u = kron3(e0 , I  , I );
  BCxy1u = kron3(e1 , I  , I );
  BCxy0v = kron3(e0 , Ih , I ); 
  BCxy1v = kron3(e1 , Il , I ); 
  BCxy0w = kron3(e0 , I  , I ); 
  BCxy1w = kron3(e1 , I  , I );

  BCzx0u = kron3( I , e0 , I ); 
  BCzx1u = kron3( I , e1 , I ); 
  BCzx0v = kron3( I , e0 , I ); 
  BCzx1v = kron3( I , e1 , I ); 
  BCzx0w = kron3( Ih, e0 , I ); 
  BCzx1w = kron3( Il, e1 , I );


  drfxy0u =  BCxy0u*xy0u;  
  drfxy1u = -BCxy1u*xy1u;   
  drfxy0v =  BCxy0v*xy0v;  
  drfxy1v = -BCxy1v*xy0v;
  drfxy0w =  BCxy0w*xy0w;   
  drfxy1w = -BCxy1w*xy0w;

  drfzx0u =  BCzx0u*zx0u;  
  drfzx1u = -BCzx1u*zx1u; 
  drfzx0v =  BCzx0v*zx0v;  
  drfzx1v = -BCzx1v*zx0v;
  drfzx0w =  BCzx0w*zx0w;   
  drfzx1w = -BCzx1w*zx0w;

  drf =[drfxy0u + drfxy1u + drfxy0v + drfxy1v + drfxy0w + drfxy1w+...
        drfzx0u + drfzx1u + drfzx0v + drfzx1v + drfzx0w + drfzx1w];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now comes the scaling

Lr = Lrf*dx^2;
Dr = Drf*dx;
Gr = Grf*dx;
Pr = Prf;
lr = Lrf*dx^2;
dr = drf*dx;








