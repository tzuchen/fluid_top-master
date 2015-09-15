
clear mixnorm;
iter =50;
tic,q20000 = PetscMaprefine(20000,iter),toc
save /home/tzuchen/Desktop/T'emp computation result'/q20000
