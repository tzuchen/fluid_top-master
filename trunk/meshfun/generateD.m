function Dr = generateD(cdinfo)
% given multi dimensional coordinate information
% this function generates the Divergent operator
%
%  Tzu-Chen Liang 10/31/2005

t0   = clock;
dim  = cdinfo.dim;

for i = 1:dim
     n     = cdinfo.coord{i}.n; 
     gsize = cdinfo.coord{i}.gsize;
     p     = cdinfo.coord{i}.p;

     I     = speye(n,n);
     E     = sparse(2:n,1:n-1,1,n,n);

     D        = I-E';
     if p == 1
         D(n,1)  = -1;
     else
         D       = D(:,2:n);
     end

     for j = 1:dim
          if i==j      
             	Op{i,j} = D/gsize;
	  else 
	  	Op{i,j} = I/gsize;          
          end

     end
end

for j= 1:dim
   Dt{j} = kronn(Op{dim:-1:1,j});
end
   Dr = cell2mat(Dt);

if cdinfo.showtime==1
  disp(sprintf('Matlab: Time to generate D : %f sec ',etime(clock,t0))); 
end
