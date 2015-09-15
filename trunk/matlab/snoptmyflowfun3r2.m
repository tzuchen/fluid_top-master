 function [F,G] = snoptmyflowfun3r2(x)
% Computes F and dense Jacobian G for the 3d flow problem.
global meshdata
global info
% x is now alpha

np        = meshdata.np;
ngd       = meshdata.ngd;
alphalist = meshdata.alphalist;


M = (info.alphamax-info.alphamin);



meshdata.betavec  = x;
meshdata          = beta2alpha2(meshdata);
meshdata.alphavec = meshdata.alphavec*M;
% Reset objective vector c
meshdata  =  ObjectiveSetting(meshdata,info);
 t  = cputime;
meshdata = diffsolverpetsc3d(meshdata,1,info);
 t  = cputime - t;
 
%%%%%%%%%%%%
if t >0.5
    backup.alphavec = meshdata.alphavec;
    backup.betavec  = meshdata.betavec;
    backup.dfdalpha = meshdata.dfdalpha;
    save('backup','backup');
    disp(sprintf('A backup alphavec and betavec is saved! (with dfdalpha)'))
   
end
%%%%%%%%%%%%
 
    
if info.alphacon == 0
   % K = 1e15;  % to make the G close to 1  (15)
    F = meshdata.obj;
    G = ((meshdata.dfdalpha)'*meshdata.M)';

    if ~isfield(meshdata,'Kscale')
        meshdata.Kscale = 1/max(max(G));
    end
    K = meshdata.Kscale;
    F = F*K;
    G = G*K;
    
    
    
else
    K = 1e-1;  % to make the G close to 1
    F = [ meshdata.obj*K; sum(x)];
    G = [ meshdata.dfdalpha(alphalist)*K*M ];
end
Obj = 1;



J = [   ];

% iGfun = J(:,1); jGvar = J(:,2); G = J(:,3);
% G = sparse(iGfun,jGvar,G);
% G = full(G);
