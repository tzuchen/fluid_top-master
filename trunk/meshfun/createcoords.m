function coordinfo = createcoords(varargin)
% Given one dimensional coordinate information, this function
% creates the multiple dimentional coordinate infomation. 
%
% Tzu-Chen Liang   10/30/2005

coordinfo.dim      = nargin; 
coordinfo.coord    = varargin;
coordinfo.showtime = 0;

fieldinfo        = cell2mat(coordinfo.coord);
coordinfo.grid   = cell2mat({fieldinfo.n});


for i = 1:coordinfo.dim
   gridspam= coordinfo.grid; 
   if coordinfo.coord{i}.p ~= 1 
        gridspam(i) = gridspam(i) -1;
   end
   nofgrid{i} = prod(gridspam);
   coordinfo.gridspam{i} = gridspam;   
end

coordinfo.nofgrid = nofgrid; 

% for pressure
coordinfo.gridspam{coordinfo.dim+1} = coordinfo.grid;
coordinfo.nofgrid{coordinfo.dim+1}  = prod(coordinfo.grid);

