function  fr  = generatefr(bodyf,cdinfo);
% given multi dimensional coordinate information
% and body force, the function generates body force.
%
% Tzu-Chen Liang 10/31/2005

t0   = clock;
dim  = cdinfo.dim;

fr = [];
for i = 1:dim
   fr  = [fr; bodyf{i}*ones(cdinfo.nofgrid{i},1)];
end

if cdinfo.showtime==1
  disp(sprintf('Matlab: Time to generate fr : %f sec ',etime(clock,t0))); 
end
