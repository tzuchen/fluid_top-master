function [S] = SparseLoad(irfile,jcfile,prfile,m,n)
%
% This function load the sparse matrix saved by Petsc
% in row compressed form (ir,jc,pr). In Petsc side, use
% the function "SparseSave"
%
% Tzu-Chen Liang 9-11-2007

      ir  = load(irfile);
      jc  = load(jcfile);
      pr  = load(prfile);

   zind = find(ir==0);

   for i = length(zind):-1:2  
      ir(zind(i):end) = ir(zind(i):end) + ir(zind(i)-1); 
   end
   ir = ir(setdiff(1:length(ir),zind(2:end)));
  
   si = zeros(1,ir(end));
   for i = 1:length(ir)-1
     si(ir(i)+1:ir(i+1)) = i;
   end
  
   ir = si;  
   jc = jc+1;

 
 if nargin >3  
   S = sparse(ir,jc,pr,m,n);
 else
   S = sparse(ir,jc,pr);
 end
