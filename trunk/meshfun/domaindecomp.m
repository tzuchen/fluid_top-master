function [p,r] = domaindecomp(sects,cdinfo)
% This function decomposes the domain into equal-sized
% sub-domains and return the permutation vector. 
%  
%   pind{ (u,v,w,p) , (x,y,z) , (sec1 sec2,..)  }
%   qind{ (u,v,w,p) , number of group }

t0       = clock;
dim      = cdinfo.dim;
nofgrid  = cdinfo.nofgrid;
gridspam = cdinfo.gridspam;

pind = [];

for i = 1:dim+1 % include pressure
     for j = 1:dim
          x  = ind2ijk(1:nofgrid{i},i,cdinfo);
          searchlist = 1:nofgrid{i};
          for k = 1:sects(j)-1   
              pind{i,j,k} = searchlist(find(x(j,searchlist)<=ceil(k*gridspam{i}(j)/sects(j))));             
              searchlist  = setdiff(searchlist,pind{i,j,k});                
          end
              pind{i,j,sects(j)} = searchlist;
     end
end

% This is bad, can only work for 3-d case!
for i = 1:dim+1
        w = 1;
        for k1 = 1:sects(1)
               for k2 = 1:sects(2)
                   for k3 = 1:sects(3)
          qind{i,w} = intersect(intersect(pind{i,1,k1},pind{i,2,k2}),pind{i,3,k3});
          w = w+1;
                   end
          
		end
	end
end    

ind = 0;
for i = 2:dim+1
       ind = ind + nofgrid{i-1};
       for w = 1:prod(sects)
	     qind{i,w} = qind{i,w}+ ind; 
       end
end


p = [qind{:,:}];
r(p) = 1:length(p);

if cdinfo.showtime==1
  disp(sprintf('Matlab: Time to decompose the domain : %f sec ',etime(clock,t0))); 
end












