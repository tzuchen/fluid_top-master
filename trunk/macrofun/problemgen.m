function [Ap,Ac,Pp,Pr,bp,p2Amatp,p,r,cdinfo] = problemgen(nproc,den,bodyf,mu)
%
% Use this function to generate problem data quickly
%
% Tzu-Chen Liang 2-25-2006


X = createcoord(0,5,20*den,1);
Y = createcoord(0,1,4*den,0);
Z = createcoord(0,1,4*den,0);

cdinfo = createcoords(X,Y,Z);
cdinfo.showtime  = 1;

Lr     = generateL(cdinfo);
Dr     = generateD(cdinfo);
Gr     = Dr';
Pr     = abs(Dr')/2/max(max(Dr));  % Be careful when large size!

% Now comes the BC setting

  BCs{1,1}  = [];               % 0yzu  
  BCs{1,2}  = [];               % 1yzu
  BCs{1,3}  = zeros(10,10);     % x0zu
  BCs{1,4}  = zeros(10,10);     % x1zu
  BCs{1,5}  = zeros(10,10);     % xy0u
  BCs{1,6}  = zeros(10,10);     % xy1u

  BCs{2,1}  = [];               % 0yzv  
  BCs{2,2}  = [];               % 1yzv
  BCs{2,3}  = zeros(10,10);     % x0zv
  BCs{2,4}  = zeros(10,10);     % x1zv
  BCs{2,5}  = zeros(10,10);     % xy0v
  BCs{2,6}  = zeros(10,10);     % xy1v

  BCs{3,1}  = [];               % 0yzw  
  BCs{3,2}  = [];               % 1yzw
  BCs{3,3}  = zeros(10,10);     % x0zw
  BCs{3,4}  = zeros(10,10);     % x1zw
  BCs{3,5}  = zeros(10,10);     % xy0w
  BCs{3,6}  = zeros(10,10);     % xy1w


[lr,dr] = generatelrdr(BCs,cdinfo);
fr      = generatefr(bodyf,cdinfo);
[A,b]   = generateAb(mu,Lr,Dr,lr,dr,fr,cdinfo);
[p,r]   = domaindecomp([2 2 2],cdinfo);
Ap      = A(p,p);
bp      = b(p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P      = [ Pr ; sparse(cdinfo.nofgrid{cdinfo.dim+1},size(Pr,2))];
Pp     = P(p,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set a preconditioner matrix Ac
% Note: now it is just an identity. When it is not just 
% an identity, DON'T forget the permutation

Ac      = speye(length(p)); 


p2vmat  = Pr;
p2Amat  = cat(1,p2vmat,sparse(size(A,1)-size(p2vmat,1),size(p2vmat,2)));
p2Amatp = p2Amat(p,:);


  
