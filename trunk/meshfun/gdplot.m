function gdplot(i,dim,cdinfo,param)
% given a set of points, plot their positions
% This function is used for debugging
%
% Tzu-Chen Liang 11/6/2005

[x,ind]  = ind2loc(i,dim,cdinfo);
if ~isempty(x)

 plot3(x(1,:),x(2,:),x(3,:),param);

end
axis([cdinfo.coord{1}.l cdinfo.coord{1}.r cdinfo.coord{2}.l cdinfo.coord{2}.r ,cdinfo.coord{3}.l cdinfo.coord{3}.r])
axis equal
xlabel('x')
ylabel('y')
zlabel('z')






