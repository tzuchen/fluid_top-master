%
%
clear M
%load Mc
%load Md
%Ma=Mc;
%Mb=Md;
load Ma
load Mb

ilist =1:300%[1 30 50 100]

for i = 1:length(ilist)
 k = ilist(i);

[Xa,mapa ] = frame2im(Ma(k));
[Xb,mapb ] = frame2im(Mb(k));
 B         = 255*ones(size(Xa,1),40,3);

imagesc([Xa B Xb])



%axis([1  size(Xa,2)+size(Xb,2) 1 size(Xa,1)])
axis equal
axis tight
%axis([1  size(Xa,2)+size(Xb,2) 1 size(Xa,1)])
M(i) = getframe

end
