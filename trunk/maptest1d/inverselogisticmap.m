function x = inverselogisticmap(y)
 
 c = 4;

 x = [0.5*(1-sqrt(1-y)), 0.5*(1+sqrt(1-y))];
