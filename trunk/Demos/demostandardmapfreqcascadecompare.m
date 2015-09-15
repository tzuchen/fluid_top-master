% demostandardmapfreqcascadecompare
% 
%
% Tzu-Chen 10-17-2007

close all

 

for i=1:20
  fname = ['wnlog_20000_b', num2str(i)];
  load(fname)
  fname = ['wnlog_20000_', num2str(i)];
  load(fname)
  fname = ['wnlog_20000_', num2str(i),'s'];
  load(fname)
end


noffile = 10

for i=[1 4 8]
  fname = ['wnlog_20000_', num2str(i)]
  w = eval(fname);
  loglog(w(:,1),w(:,2),'r')
  hold on

  fname = ['wnlog_20000_', num2str(i),'s']
  w = eval(fname);
  loglog(w(:,1),w(:,2),'g')
  hold on

  fname = ['wnlog_20000_b', num2str(i)]
  w = eval(fname);
  loglog(w(:,1),w(:,2),'b')
  hold on

end
