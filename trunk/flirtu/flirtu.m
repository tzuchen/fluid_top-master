function [X,y,S,option] = flirtu(C,A,b,option)
% An implementation of Primal-Dual SDP solver
% It solves the SDP problem in the standard form:   
%         min  trace(C*X)
%         s.t. trace(A{i}*X) = b(i)
%              X is PSD
% Where A{i} is a n*n cell matrix with length m
%       C    is a n*n matrix
%       b    is a vector with length m
%
% Option list      
%    option.blkstruc : The block structure of X. It is a vector
%                     
%    option.sol : an initial feasible solution, if this is not 
%                 given, Big-M method is applied to generate 
%                 one. It is given by 
%                           sol.X0 = X0
%                           sol.y0 = y0  
%                           sol.S0 = S0   
%     
%    option.etol : the algorithm stops when the duality gap is less than it
%    option.alpha: alpha parameter for line search
%    option.beta : beta parameter for line search
%
                                      
%    
%
% Tzu-Chen Liang           12/3/2005

if nargin == 3
   option.m  = length(b);
   option.n  = size(C,1); 
   [X,y,S,option]   = flirtu(C,A,b,option);
else
   if ~isfield(option,'m')
       option.m  = length(b);       
   end
   if ~isfield(option,'n')
       option.n  = size(C,1);
   end
   if ~isfield(option,'blkstruc')
       option.blkstruc = [size(C,1)];     
   end

   if ~isfield(option,'sol')
      % No initial point given
       option.blkstruc = [option.blkstruc 1 1 ]; 
       X0              = speye(option.n);
       y0              = zeros(option.m,1);
       S0              = speye(option.n);
       option.m        = option.m + 1;
       option.n        = option.n + 2;

       M               = 10^2;
       N               = 10^2;     
       Cm              = blkdiag(C,M,0);
       for i = 1:option.m-1
            Am{i} = blkdiag(A{i},b(i)-trace(A{i}*X0),0);
       end
            Am{option.m} = blkdiag(-S0+C,0,-1);
            bm           = [ b ; -N ];
       u                 = N - trace((S0-C)*X0);
       v                 = M;
  
       option.bigM      = 1; 
       option.sol.X0    = blkdiag(X0,1,u);
       option.sol.y0    = [y0; 1];
       option.sol.S0    = blkdiag(S0,v,1);    

       option.A         = Am;
       option.C         = Cm;
       option.b         = bm;
 

   end
   
   if ~isfield(option,'beta')
       option.beta = 0.9;
   end
   if ~isfield(option,'alpha')
       option.alpha = 0.5;
   end
   if ~isfield(option,'etol')
       option.etol  = 1e-6;
   end

   if ~isfield(option,'parallel')
       option.parallel = 1;
   end

end

if and(nargin==4,isfield(option,'sol'))

     if isfield(option,'A')
          A = option.A;
          C = option.C;
          b = option.b;
     end

     X = option.sol.X0;
     S = option.sol.S0;
     y = option.sol.y0;

     m        = option.m;
     n        = option.n;
     blkstruc = option.blkstruc;    
     gama     = n/(n+sqrt(n));
     [index]  = blkindex(blkstruc);
     nv       = length(index);

     Amat     = sparse(m,nv);
     for i = 1:m
         Amat(i,:) = A{i}(index);
     end

     Inv      = speye(nv);
     Im       = speye(m);
     Onv      = sparse(nv,nv);
     Om       = sparse(m,m);
     Omnv     = sparse(m,nv);
     
     dgap     = trace(X*S);
     count    = 0;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if option.parallel ==1
      nproc   = 5;
      petscopt.ksp_type = 'cg';
      PetscInitialize(nproc,'sdplinkernal',petscopt);
      socketp = openport();
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


     while dgap > option.etol
        count = count + 1;
     	Xh  = sqrtm(full(X));
     	D   = Xh*sqrtm(inv(Xh*S*Xh))*Xh;
     	Dh  = sqrtm(D); 
     	mu  = trace(X*S)/n;
     	R   = gama*mu*inv(X)-S;
     	Rp  = Dh*R*Dh;

     	Apmat = sparse(m,nv);  
     	for i = 1:m
              temp = Dh*A{i}*Dh;
              Apmat(i,:) = temp(index);
        end

        Nmat  = [Onv   Apmat' Inv  ;
                 Apmat Om     Omnv ;
                 Inv   Omnv'  Inv  ];
        Rhs   = [zeros(nv,1); zeros(m,1); Rp(index) ];

   Imat = speye(size(Nmat,1));
 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if option.parallel == 1 
      send(socketp,count);   
      send(socketp,Nmat);
      send(socketp,Imat);
      send(socketp,Rhs);
      dxys = receive(socketp);
   else
      dxys  = Nmat\Rhs; % Use matlab solver
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
        dxp   = dxys(1:nv);
        dy    = dxys(nv+1:nv+m); 
        dsp   = dxys(nv+m+1:end);

     	dXp = sparse(n,n);
     	dXp(index) = dxp;
     	dXp = dXp+dXp';
     	dX  = Dh*dXp*Dh;
     
        ds    = -Amat'*dy;
	dS  = sparse(n,n);
        dS(index) = ds;
        dS = (dS + dS'- diag(diag(dS)));

        t = bklinesearch(X,S,dX,dS,option.alpha,option.beta,n,sqrt(n));

        X = X+t*dX;
        S = S+t*dS;
        y = y+t*dy; 



        dgapn = trace(C*X)- b'*y;
        dgap = [dgap  dgapn];
        disp(sprintf('iter = %i, duality gap = %d', count ,dgapn))

   end

   if option.bigM == 1;
       X = X(1:end-2,1:end-2);
       S = S(1:end-2,1:end-2);
       y = y(1:end-1);
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if option.parallel ==1  
       PetscFinalize(socketp);
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end





