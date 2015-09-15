

load Mb
figure 
hold on

for i = 0:9

   X = Mb(10*i+1).cdata;
   z = 10*i 
  
   surface('XData',[0 1;0 1],...
      'YData',[0 0;1 1],...
      'ZData',[z z; z z],...
      'CData',X,...
      'FaceColor','texturemap');

end
