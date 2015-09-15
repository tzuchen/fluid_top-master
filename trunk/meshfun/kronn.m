function M = kronn(varargin)


M = varargin{end};
for i = nargin:-1:2
    M = kron(varargin{i-1},M); 
end
