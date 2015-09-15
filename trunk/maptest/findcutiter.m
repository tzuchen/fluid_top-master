function [cutiter]=findcutiter(varn,cutsigma)

[nofcase,nofiter] = size(varn);

for i =1:nofcase
   cutiter(i) = interp1(varn(i,:),0:nofiter-1,cutsigma);
end
