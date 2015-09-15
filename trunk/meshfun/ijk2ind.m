function ind = ijk2ind(ijk,dir,cdinfo)
% Given the ijk coordinates and direction, this function
% returns the one dimensional location.
%
%  Tzu-Chen Liang 11/15/2005

gridspam = cdinfo.gridspam{dir};
dim      = cdinfo.dim;
ind      = ijk(1,:);

for i = 2:dim
   ind = ind + (ijk(i,:)-1)*prod(gridspam(1:i-1)); 
end
