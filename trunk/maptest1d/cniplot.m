function cniplot(n,l)

clist= zeros(l,1);

clist(1) = 1;

for i = 2:l
 clist(i) = 1/n+clist(i-1)*(n-2)/n;
end

llist = 1:l;
%llist = log(llist);
llist = llist./(0.25*n*log(n));

%plot(llist,clist)
semilogx(llist,clist)
%plot(llist,clist)


