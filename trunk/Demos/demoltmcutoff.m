% demoltmcutoff

n = 1500;
nit = 500;

varp=[];
[A2] = maprefine2(n,[],@twistmap2,0.4,0.4,0.3,1);
[A1] = maprefine2(n,[],@twistmap2,0.6,0.6,0.3,1);




difflist = [0 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 ];

for i = 1:length(difflist)

 [Mdiff]  = fourierdiffuseM(n,difflist(i),1,1);
 [Xa1,varx,vary,h0norm] = testltm(A1.A,A2.A,Mdiff,'cosx',nit,0);
 varp(i,:) = varx;
 save /home/tzuchen/Desktop/T'emp computation result'/demoltmcutoffdata 

end




load /home/tzuchen/Desktop/T'emp computation result'/demoltmcutoffdata 

close all
figure
hold on

varp = varp([1, 3:end],:);
h = plot(varp')
set(h,'linewidth',2)
grid on
box on
xlabel('iterations')
ylabel('...')
legend('D=D^*','D=1e-7','D=1e-6','D=1e-5','D=1e-4','D=1e-3')


