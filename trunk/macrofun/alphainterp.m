function [al,w2]= alphainterp(fname,cdinfo,nx,ny,nz)
% This routine do the interpolation of alpha*files
% and double the resolution
% 
%  Tzu-Chen Liang 

ns = cdinfo.gridspam{4};
nxf = ns(1);
nyf = ns(2);
nzf = ns(3);


load(fname);

w = reshape(eval(fname),nx,ny,nz);

xr = cdinfo.coord{1}.r;
yr = cdinfo.coord{2}.r;
zr = cdinfo.coord{3}.r;

dxf = (cdinfo.coord{1}.r-cdinfo.coord{1}.l)/nxf;
dx  = (cdinfo.coord{1}.r-cdinfo.coord{1}.l)/nx;

% make it periodical
wp = w([end,1:end,1],:,:);


[XI,YI,ZI] = meshgrid(dxf/2:dxf:yr,dxf/2:dxf:xr,dxf/2:dxf:zr);
[X,Y,Z] = meshgrid(dx/2:dx:yr,-dx/2:dx:xr+dx/2,dx/2:dx:zr);


w2 = interp3(X,Y,Z,wp,XI,YI,ZI,'linear',0);

al = reshape(w2,nxf*nyf*nzf,1);

send2file('beinpfiled',al);

disp(sprintf('File %s has size %d %d %d is loaded', fname,nx,ny,nz ))
disp(sprintf('The above file is interpolated and has new size %d %d %d',nxf,nyf,nzf ))
disp(sprintf('and saved as %s ', 'beinpfiled' ))


