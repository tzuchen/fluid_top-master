 function [x,F,inform] = snoptflow3r2(meshdata,x)
% First derivatives are provided.
 global meshdata
 global info

snoptflow.spc = which('snoptflow.spc');
snprintfile('snoptflow.out');
snsummary  ('snoptflow.sum');
snspec     ('snoptflow.spc');
snseti     ('Major Iteration limit', 1250);

%Get condensed data.
[x0,xlow,xupp,Flow,Fupp,A,iAfun,jAvar,iGfun,jGvar] = myflowr2();

% Max
 snset('Maximize');
%snset('Minimize');

[x,F,inform] = snopt(x,xlow,xupp,Flow,Fupp,'snoptmyflowfun3r2', ...
	    	         A, iAfun, jAvar, iGfun, jGvar);

                 
%  meshdata.betavec = x;
%  meshdata = beta2alpha(meshdata);
%  x   =  meshdata.alphavec;                    
%  x = x *(info.alphamax - info.alphamin);                 
                     
snset('Defaults');
snsummary off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function [x,xlow,xupp,Flow,Fupp,A,iAfun,jAvar,iGfun,jGvar] = myflowr2()

 global meshdata
 global info

n      =  meshdata.nbeta;   % number of variables

if info.alphacon == 0 
 
     neF           =  1; 
     [iGfun,jGvar] = find(ones(neF,n));
     iGfun         = iGfun';
     jGvar         = jGvar';
     iAfun         = [];  
     jAvar         = [];
     A             = [];
     
     Flow(1)       = -Inf; 
     Fupp(1)       =  Inf;
else
    
    neF            = 2;
    totalmat       = info.matratio*(1 - 0)*meshdata.np;
    [iGfun,jGvar]  = find(ones(1,n));
    iGfun          = iGfun';
    jGvar          = jGvar';    
    iAfun          = [2*ones(1,n)]';  
    jAvar          = [1:n]';  
    A              = [ ones(1,n)'];

    Flow           = zeros(neF,1);
    Fupp           = zeros(neF,1);

    Flow(1)        = -Inf; 
    Fupp(1)        =  Inf;
    Flow(2)        =  totalmat;    
    Fupp(2)        =  inf;%totalmat; 
end


ObjAdd  = 0;

% Ranges for x.
xlow = 0*ones(n,1);  
xupp = 1*ones(n,1);


x      = zeros(n,1);
xstate = zeros(n,1);
xmul   = zeros(n,1);

Fmul   = zeros(neF,1);
Fstate = zeros(neF,1);

