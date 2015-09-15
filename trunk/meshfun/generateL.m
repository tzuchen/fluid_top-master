function Lr = generateL(cdinfo)
% given multi dimensional coordinate information
% this function generates the Laplacian operator
%
% Tzu-Chen Liang 10/31/2005

t0       = clock;
dim      = cdinfo.dim;
gridspam = cdinfo.gridspam;
rpflag   = 0;

for i = 1:dim
     n     = cdinfo.coord{i}.n; 
     gsize = cdinfo.coord{i}.gsize;
     p     = cdinfo.coord{i}.p;

     I     = speye(n,n)/gsize;
     e  	 = ones(n,1);
     L     = spdiags([e -2*e e], -1:1, n, n)/gsize;
     if p==1
     	  L(1,n) = 1/gsize;
   	  L(n,1) = 1/gsize;  
     end
 
     for  j = 1:dim
         if i==j
           Op{i,j} = L;
         else
           Op{i,j} = I;
         end
     end  
end

for i = 1:dim
   Lt{i} = sparse(prod(gridspam{i}),prod(gridspam{i}));
end


for j = 1:dim
    for i = 1:dim 
    	if cdinfo.coord{i}.p == 0  
            n     = cdinfo.coord{i}.n; 
   	    rpflag  = 1; 
    	    tempOp  = Op{i,j};
    	    Op{i,j} = Op{i,j}(1:n-1,1:n-1); 
    	end

    	Lt{i} = Lt{i} + kronn(Op{dim:-1:1,j});
 	if rpflag == 1;
    		Op{i,j} = tempOp;
        	rpflag = 0;
   	end
    end
end

Lr =blkdiag(Lt{:});

if cdinfo.showtime==1
  disp(sprintf('Matlab: Time to generate L : %f sec ',etime(clock,t0))); 
end

