clear all
%
% fd = 1: solve the problem
% fd = 0: save the problem data
%
% plotstructure = 1: plot it
% plotstructure = 0: don't plot
%
fd = 0;
plotstructure = 0; 

den  = 20;
dx = 1/den;
dy = 1/den;

X = createcoord(0,1,den,0);
Y = createcoord(0,1,den,0);
nx = (den-1+X.p)*den;
ny = den*(den-1+Y.p);

cdinfo = createcoords(X,Y);
cdinfo.showtime  = 1;

Lr     = generateL(cdinfo);
Dr     = generateD(cdinfo);
Gr     = Dr';
Pr     = abs(Dr')/2/max(max(Dr));  % Be careful when large size!



% Now comes the BC setting
lr = zeros(size(Lr,1),1);
dr = zeros(size(Dr,1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ygu,xgu] = meshgrid(dx:dx:1-dx,dx/2:dx:1-dx/2); 
[ygv,xgv] = meshgrid(dx/2:dx:1-dx/2,dx:dx:1-dx);

xbc1 = zeros(size(Lr,1),1);
xbc2 = zeros(size(Lr,1),1);


lb = 0.2;
ub = 0.8;
ind1 = find(and(and(xgu<=ub,xgu>=lb),ygu<=dy));
xbc1(ind1) = 1*(ub-xgu(ind1)).*(lb-xgu(ind1));
ind2 = nx+find(and(and(xgv<=ub,xgv>=lb),ygv<=dy));
xbc1(ind2) = 0;

lb = 0.2;
ub = 0.8;

ind3 = find(and(and(xgu>=lb,xgu<=ub),ygu>=1-2*dy));
xbc2(ind3) = 1*(ub-xgu(ind3)).*(lb-xgu(ind3));
ind4 = nx+find(and(and(xgv>=lb,xgv<=ub),ygv>=1-2*dy));
xbc2(ind4) = 0;

xbc = xbc1/sum(xbc1) + xbc2/sum(xbc2);

ind = [ind1; ind2; ind3; ind4];

lr  = Lr*xbc;
dr  = Dr*xbc;
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bodyf{1}    =  0;
bodyf{2}    =  0;


mu = 1;

fr      = generatefr(bodyf,cdinfo);

solind = setdiff(1:size(Lr,1),ind);
Lrs  = Lr(solind,solind);
Drs  = Dr(:,solind);
drs  = dr;
lrs  = lr(solind);
frs  = fr(solind);
[Ar,br]   = generateAb(mu,Lrs,Drs,lrs,drs,frs,cdinfo);

[A,b]   = generateAb(mu,Lr,Dr,lr,dr,fr,cdinfo);


sol = Ar\br;
x(solind) = sol(1:length(solind));
x(ind)=xbc(ind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%testmeshfun
n = size(Lrs,1);
m = size(Drs,2);
%br = b(1:n)';
M  = abs(Drs);

%cvx_begin
%    variable x(n)
%    minimize( x'* (-Lrs) * x -lrs'*x )
%    subject to
%        Drs * x == drs;
       % norm(x,2) <= 0.01;
%cvx_end
%sol = x;
%clear x
%x(solind) = sol(1:length(solind));
%x(ind)=xbc(ind);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



vx = x(1:nx);
vy = x(nx+1:end);

Vx = reshape(vx,den-1+X.p,den);
Vy = reshape(vy,den,den-1+Y.p);

if X.p ==0
Vx = [zeros(den,1)'; Vx; zeros(den,1)'];
else
Vx = [ Vx; Vx(1,:)];
end

Vy = [zeros(den,1)   Vy zeros(den,1)];

Vx = 0.5*(Vx(1:end-1,:)+Vx(2:end,:));
Vy = 0.5*(Vy(:,1:end-1)+Vy(:,2:end));


[gX,gY] = meshgrid(dx/2:dx:1-dx/2,dx/2:dx:1-dx/2);
Ck =2*1/( max(abs(x))*den);

close 
figure
hold on
for i = 1:(den)^2
 
  line([gX(i),gX(i)+Ck*Vx(i)],[gY(i),gY(i)+Ck*Vy(i)])
 % line([gY(i),gY(i)+Ck*Vy(i)],[gX(i),gX(i)+Ck*Vx(i)])
end











