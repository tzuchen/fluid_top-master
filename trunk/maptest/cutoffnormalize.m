function cutoffnormalize(var2,noval)



[n, iter] = size(var2);

for i = 1:n
  for j = 1:iter-1
    if var2(i,j+1)>=var2(i,j);
      var2(i,j+1) = var2(i,j)-1e-9; 
    end
  end
end


for i = 1:n
  itern(i)   = interp1(var2(i,:),0:iter-1,noval);
end




%colorlist = 'gcmybgcmybgcmyb'
figure
hold on

for i = 1:n

h(i) = plot([0:iter-1]/itern(i),var2(i,:));
  
end
axis([0 3 0 1])
set(h,'linewidth',2)
grid on
box on
