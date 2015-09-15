function [index] = blkindex(blkstruc)
% Given the block structure of a symmetric matrix,
% this function generate the index to the significant  
% variables

 X = [];
 for i = 1:length(blkstruc)
    X = blkdiag(X,sparse(triu(ones(blkstruc(i)))));
 end

index  = find(X);
