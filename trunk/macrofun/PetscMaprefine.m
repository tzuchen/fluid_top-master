function [mixnorm1,mixnorm2,X] = PetscMaprefine(n,iter,doSmoothing,doFFT,FFTdiffu,IC,num)
%
%  Tzu-Chen Liang 9-29-2006

nproc = 3;

option.withMatlab = 1;


if nargin>2
  if doSmoothing >0
    option.doSmoothing = doSmoothing;    
  end
end
if nargin>3
  if doFFT==1
    option.doFFT = 1;
  else 
    option.doFFT = 0; 
  end
end
if nargin>4
  if FFTdiffu>=0
     option.doFFTdiffusion= FFTdiffu;
  end
end

if nargin>5
     option.IC= IC;
end

if nargin>6
  MaprefineVer = ['Maprefine',num2str(num)];
else
  MaprefineVer = ['Maprefine'];
end




PetscInitialize(nproc,MaprefineVer,option);
socketp = openport();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEND part %
 send(socketp,n);
 send(socketp,iter);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RECEIVE part %
x0      = receive(socketp);
mixnorm1 =  receive(socketp);
mixnorm2 =  receive(socketp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mixnorm1 = mixnorm1(1:iter);
mixnorm2 = mixnorm2(1:iter);

PetscFinalize(socketp);


X1 = reshape(x0,n/2,n);
X2 = X1(end:-1:1, end:-1:1);
X  = [X1;X2];
imagesc(X);
