% Test routine to generate a Markov matrix with arbitrary diffusion.
%
%
% Tzu-Chen Liang 2/1/2006

% need variables:
%                nofline 
%
%
% 

nofline = length(y);
onevec           = ones(nofline,1);
%onevec(find(y0)) = 0;
onevec           = onevec(1:end-1);

M = diag(sparse(onevec),-1);
d = sparse(nofline,nofline);
G = sparse(nofline,nofline);

for i = 1: nofline
    for j = i: nofline
         dij = norm([y(i)-y(j),z(i)-z(j)]);
         if and(dij<0.03,0<dij)
          d(i,j) = dij;
          d(j,i) = dij;
          G(i,j) = exp(-(dij*60).^2)/pi;
          G(j,i) = G(i,j);
         end

   
    end
    G(i,i) = 1/pi;
end



%W = G*M;
%W2 = W+M;
%for i = 1:nofline
%  W2(i,:) = W2(i,:)/sum(W2(i,:));
%  G(i,:)  = G(i,:)/sum(G(i,:));
%end

%reso = sort(eig(full(W2)))
