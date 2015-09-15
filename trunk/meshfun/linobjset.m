function c = linobjset(cpos,cdinfo);
% This function sets the linear objective function
% it returns a cell vector c, which is just df/dx
% To use it as a rhs vector just do 
%  
%           c = cell2mat(c')
%
% Tzu-Chen Liang  11/8/2005

t0  = clock;
n   = size(cpos,2);
dim = cdinfo.dim;

for i = 1:dim+1
 ct{i} = zeros(cdinfo.nofgrid{i},1);
end


for setn = 1:n
  val      = cpos{setn}(end);
  dir      = cpos{setn}(end-1);
  ngd      = cdinfo.nofgrid{dir};
  dirgrid  = ind2loc(1:ngd,dir,cdinfo);
  indexh   = 1:ngd;
      
      for i = 1:dim         
          lb = cpos{setn}(2*i-1);
          ub = cpos{setn}(2*i  );
          indexl = indexh(find(lb<dirgrid(i,indexh)));
          indexh = indexl(find(dirgrid(i,indexl)<ub));
      end 
      cind{setn} = indexh;    
      ct{dir}(cind{setn}) = val; 
end

c = ct;
%c = cell2mat(ct');

if cdinfo.showtime==1
  disp(sprintf('Matlab: Time to set c : %f sec ',etime(clock,t0))); 
end
