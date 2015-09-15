%democutoff
%
%
% 
% Plot a sequence of trajectories that show the typical cutoff
%
% Tzu-Chen Liang 6-19-2007



klist = 0:1:200
nlist = 0:10:100
clear x
i = 1;
for n = nlist
   
  x(i,:) = exp(-0.1*exp(0.1*(klist-n)));
  i=i+1;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
figure 
hold on

 for i=1:size(x,1)
     h=plot(x(i,:))
     set(h,'linewidth',2)
 end

axis([0 150 0 1])
grid on
box on
xlabel('iterations')
ylabel('measure of certainty')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cutoffnormalize(x,0.5)


grid on
box on
xlabel('normalized iterations')
ylabel('measure of certainty')
