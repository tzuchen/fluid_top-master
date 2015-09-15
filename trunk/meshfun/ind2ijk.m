function ind  = ind2ijk(i,dim,cdinfo)
% Given and nD coordinate system cdinfo 
% and an index i in the certain dim,
% return the i,j,k location in space.
% Notice : this function can work in 
% n-dimensional coordinate system.  
%    
%  Tzu-Chen Liang  11/1/2005

gridspam = cdinfo.gridspam{dim};
if size(i,2)<size(i,1) 
   i = i'; % make sure i is a row vector
end

for ndim = 1:length(gridspam)-1
    ind(:,ndim) = mod(i,gridspam(ndim));
    i       = floor((i-1)./gridspam(ndim));
    i       = i + 1; 
end

ind = [ind' ;i];

for ndim = 1:length(gridspam)-1
    ind(ndim,find(ind(ndim,:)==0)) = gridspam(ndim);
end
