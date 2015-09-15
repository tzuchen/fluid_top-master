function Pd = freqmapaddD(P,n,statelist,Dstar,D);
%  Given P with effective diffusion rate Dstar
%  This function calculate Pd which has diffision D
%  note : D > Dstar   
%
%  Tzu-Chen Liang 7-28-2006


[Mdiff]  = fourierdiffuseM(n,D-Dstar,1,1);
M        = Mdiff.M(statelist);
Pd       = P*diag([M; M]); 
