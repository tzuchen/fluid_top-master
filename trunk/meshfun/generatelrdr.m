function [lr,dr] = generatelrdr(BCs,cdinfo)

t0       = clock;
dim      = cdinfo.dim;
NofBC    = 2*dim^2;
gridspam = cdinfo.grid;

for i = 1:dim
    for j = 1:dim
        gridspam      = cdinfo.grid; 
        if cdinfo.coord{i}.p~=1
	  gridspam(i) = gridspam(i) -1; 
        end
        gridspam(j)     = 1;
        BCsize{i,2*j-1} = gridspam(find(gridspam~=1));
        BCsize{i,2*j  } = BCsize{i,2*j-1};
    end
end


for i = 1:dim
    for j = 1:2*dim
       if isempty(BCs{i,j})
           xi      = BCsize{i,j}(1);
           yi      = BCsize{i,j}(2);
           BCsinter{i,j} = zeros(xi,yi);
           BCvec{i,j}    = reshape(BCsinter{i,j},xi*yi,1); 
       else
       	   [x,y]   = size(BCs{i,j}); 
           xi      = BCsize{i,j}(1);
           yi      = BCsize{i,j}(2);

           BCsinter{i,j} = interp2(BCs{i,j},1+[0:xi-1]/(xi-1)*(x-1),1+[0:yi-1]'/(yi-1)*(y-1));
           BCvec{i,j}    = reshape(BCsinter{i,j},xi*yi,1);
       end
    end
end


for i = 1:dim
    Op{i}      = speye(cdinfo.grid(i))/cdinfo.coord{i}.gsize;  % Scaling here!
    e0{i}      = zeros(cdinfo.grid(i),1); 
    e0{i}(1)   = 1;
    e1{i}      = zeros(cdinfo.grid(i),1); 
    e1{i}(end) = 1;

end

for i = 1:dim    
    for j = 1:dim     
          Opt        = Op;
          Opt{j}     = e0{j};

          if cdinfo.coord{i}.p~=1
              if i~=j
                 Opt{i} = Opt{i}(1:end-1,1:end-1);
              else
                 Opt{i} = Opt{i}(1:end-1,1);              
              end
          end
             lrt{i,2*j-1} = kronn(Opt{dim:-1:1})*BCvec{i,2*j-1}; 

          Opt        = Op;
          Opt{j}     = e1{j}; 

          if cdinfo.coord{i}.p~=1
              if i~=j
                 Opt{i} = Opt{i}(1:end-1,1:end-1);
              else
                 Opt{i} = Opt{i}(2:end,1);              
              end
          end          
             lrt{i,2*j  } = kronn(Opt{dim:-1:1})*BCvec{i,2*j  };  
   
    end
end

lr = cell2mat(lrt(:,1));
for i = 2: 2*dim 
  lr = lr + cell2mat(lrt(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dr


for i = 1:dim    
    for j = 1:dim     
          Opt        = Op;
          Opt{j}     = e0{j};
          if and(i~=j,cdinfo.coord{i}.p~=1)
             Opt{i}     = Opt{i}(:,1:end-1);         
          end
             drt{i,2*j-1} = kronn(Opt{dim:-1:1})*BCvec{i,2*j-1}; 
         
          Opt        = Op;
          Opt{j}     = e1{j};
          if and(i~=j,cdinfo.coord{i}.p~=1)
             Opt{i}     = Opt{i}(:,2:end);        
          end
             drt{i,2*j  } = -kronn(Opt{dim:-1:1})*BCvec{i,2*j  } ;

    end
end

dr = zeros(size(drt{1,1}));

for i = 1:dim
  for j = 1:2*dim;
	dr = dr + drt{i,j};
  end
end

if cdinfo.showtime==1
  disp(sprintf('Matlab: Time to generate lr and dr : %f sec ',etime(clock,t0))); 
end






