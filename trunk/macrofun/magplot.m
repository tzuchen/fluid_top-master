function wn = magplot(Xf);


Xf = abs(Xf);

n  = size(Xf,1);

[xl,yl] = meshgrid([0:n/2 n/2-1:-1:1],[0:n/2 n/2-1:-1:1]);
zl = floor((xl.^2+yl.^2).^0.5);
wn = [];

%zl = zl(1:n/2,:);
%Xf = Xf(1:n/2,:);

for i = 0:n

   ind = find(zl==i);
   if length(ind)>0
   wn(i+1) = sum(Xf(ind));
   else
   wn(i+1) = 0;
   end

end


