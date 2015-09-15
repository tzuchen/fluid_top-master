%demotentmapcutoff
%
% Tzu-Chen Liang 11-11-2007
close all
figure
hold on
box on
grid on
obj = 1

if obj == 1
 cm = jet(6);
 linestylelist= {'-','-','-','-','-','-'};
else
 cm = zeros(6,3);
  linestylelist= {'-','-','-','-','-','-'};
end


k=15
nlist  = 3:2:15;

for n=nlist 
  for i =0:k
    if i<=n-1
      mu(i+1) =  1-2^(1+i-n);
    else
      mu(i+1) = 0;
    end   
  end

   h =  plot(0:k,mu)


end
 
 
