function [ijk,rem] = loc2ijk(loc,dir,cdinfo)
%  Given the location, this function returns its i,j,k
%  coordinate and the remainder.
%    
%  Tzu-Chen Liang  11/14/2005

dim = cdinfo.dim;

for i =1:dim
    gsize = cdinfo.coord{i}.gsize;
    
    if i==dir
        temp     = loc(i,:) + gsize; 
        ijk(i,:) = floor(temp/gsize);
        rem(i,:) = temp - ijk(i,:)*gsize;
        if cdinfo.coord{i}.p ~= 1 
             ijk(i,:) = ijk(i,:) - 1;
        end
    else 
        temp     = loc(i,:)+0.5*gsize;
        ijk(i,:) = floor(temp/gsize);       
        rem(i,:) = temp - ijk(i,:)*gsize;

    end
end

