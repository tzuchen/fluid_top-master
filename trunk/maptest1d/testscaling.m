
z = 0:1/10000:1;

n = 4


%plot(wmun(z,n,n-1)-wmun(z,n,n))

PL = wmun(z,n,n-1);
PR = wmun(z,n,n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We have (PLA+PLB) = (PRA+PRB)
PLA =  wmun((1+z)/2,n-1,n-1); 
PLB =  wmun((1-z)/2,n-1,n-1);

PRA =  wmun((1+z)/2,n-1,n); 
PRB =  wmun((1-z)/2,n-1,n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We have PLA = PLA1+PLA2
PLA1 = 1./4./sqrt((1-z)/2).*( wmun( (1-sqrt((1-z)/2))/2 ,n-2,n-1) );
PLA2 = 1./4./sqrt((1-z)/2).*( wmun( (1+sqrt((1-z)/2))/2 ,n-2,n-1) );
%We have PLB = PLB1+PLB2
PLB1 = 1./4./sqrt((1+z)/2).*( wmun( (1-sqrt((1+z)/2))/2 ,n-2,n-1) );
PLB2 = 1./4./sqrt((1+z)/2).*( wmun( (1+sqrt((1+z)/2))/2 ,n-2,n-1) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We have PRA = PRA1
PRA1 = 1./4./sqrt((1-z)/2).*( wmun( (1-sqrt((1-z)/2))/2 ,n-2,n) );
% We have PRB = PRB1
PRB1 = 1./4./sqrt((1+z)/2).*( wmun( (1-sqrt((1+z)/2))/2 ,n-2,n) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We have PLA1+PLA2+PLA3+PLA4 = PRA1+PRA2;
%plot(z,PLA1+PLA2+PLB1+PLB2-PRA1-PRB1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
figure 
hold on


%plot(z, PLA1+PLA2)
%plot(z, PRA1,'r')
%plot(z, PLB1/PLB1(1)*PRB1(1),'r.')
%plot(z, PRB1) 

%semilogy(z,[PLA1;PLA2;PLB1;PLB2])
%h= semilogy(z,[PRA1;PRB1])
%set(h,'linewidth',2)

plot(z,(PLA2+PLB2)./(PLA1+PLB1))




