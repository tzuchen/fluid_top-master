function F = my2dDFT(n,freqlist,dir)
%
%  generate discrete 2-d Fourier transform matrix 
%
%  Tzu-Chen Liang 8-1-2006
%

m = length(freqlist);
freqlist = reshape(freqlist,m,1);

[r,s] = ind2sub([n,n],freqlist);

if dir >0 
 F = zeros(m,n^2);
    [p,q] = ind2sub([n,n],1:n^2);
    [r,s] = ind2sub([n,n],freqlist);
    p = p-1;
    q = q-1;
    r = r-1;
    s = s-1;

    k = -2*pi*i/n;
    F = (ones(m,1)*p).*(r*ones(1,n^2))+(ones(m,1)*q).*(s*ones(1,n^2));
    F = exp(k.*F);   

   


    %F = exp(-2*pi*i*((p-1).*(r-1)/n+(q-1).*(s-1)/n));
     

elseif dir <0
 F = zeros(n^2,m);
    [p,q] = ind2sub([n,n],1:n^2);
    [r,s] = ind2sub([n,n],freqlist);
    p = p-1;
    q = q-1;
    r = r-1;
    s = s-1; 


    k = 2*pi*i/n;
    F = (p'*ones(1,m)).*(ones(n^2,1)*r')+(q'*ones(1,m)).*(ones(n^2,1)*s');
    F = exp(k.*F); 


else
F=[];
end
 
