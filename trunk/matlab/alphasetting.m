% initial alpha setting
function meshstruc = alphasetting(meshstruc,info);

nalpha    = meshstruc.nalpha;
alphalist = meshstruc.alphalist;
np        = meshstruc.np;

matratio = info.matratio; 
alphamax = info.alphamax;
alphamin = info.alphamin;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alphar = ones(nalpha,1)* matratio*(alphamax - alphamin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meshstruc.alphavec(alphalist) = alphar;


