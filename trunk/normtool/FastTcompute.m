function [T,ilist,jlist,k] = FastTcompute(n,epsilon)

K = 1e-3;
M = [2 1; 
     1 1 ];
%K = 0.1*2*pi;
%M = [1 1;
%     0 1];



x = -n:n;
y = -n:n;

[k1,k2] = meshgrid(x,y);

k = [k1(1:end);k2(1:end)];


[ksort,ind]=sort(k(1,:).^2+k(2,:).^2);

k(1,:) = k(1,ind);
k(2,:) = k(2,ind);

s=length(ind);

kM = M'*k;
ksum = k(1,:)+k(2,:);

jlist = zeros(1,4*s); 
ilist = zeros(1,4*s);
val   = zeros(1,4*s);
nzind = 1;
q = k; 

for row = 1:s
    indq2  = find(kM(2,row)-q(2,:)==0);  % Q2 = 0 list;
    %kMs    = kM(:,indq2);
    qs     = k(:,indq2);
    Q1     = kM(1,row)-qs(1,:);
 
    indq11  = find(abs(Q1)<15);
    indq12  = find(abs(ksum(indq2)*K)-0.9*abs(Q1)>0);
    indq1   = union(indq11,indq12);

    %indq1  = 1:length(indq2);
    indrow = indq2(indq1);

    nnzrow = length(indrow);

    jlist(nzind:nzind+nnzrow-1) = indrow;
    ilist(nzind:nzind+nnzrow-1) = row*ones(1,length(indrow));
   if ~isempty(indq1)
    w = exp(-epsilon*(q(1,indrow).^2+q(2,indrow).^2));
    val(nzind:nzind+nnzrow-1) =w.* (i.^Q1(indq1)).*besselj(Q1(indq1),ksum(row)*K);
   end 
   nzind = nzind+nnzrow;
    


end
val= val(1:nzind-1);
ilist= ilist(1:nzind-1);
jlist= jlist(1:nzind-1);

T =  sparse(ilist,jlist,val,s,s);
