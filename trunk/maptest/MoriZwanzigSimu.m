function [Pbar,ErrX]=MoriZwanzigSimu(P,m,nit,x0)
% 
% this function only works for pi= uniform
%
%
% Tzu-Chen Liang  2-7-2007

n = sqrt(size(P,1)); 

B = Bgenerate(n,m);

PsiY = B;
PsiX = diag(sparse(1./sparse(sum(B))))*B';

Pbar  =PsiX*P*PsiY;

Pp  = (speye(n^2,n^2) - PsiX'*PsiY');


Errterm = reshape(x0,n^2,1);
for i = 1:nit+1
 Errterm = P'*Pp*Errterm;
end
 Errterm = PsiY'*Errterm;

ErrX = reshape(Errterm,m,m);
  
