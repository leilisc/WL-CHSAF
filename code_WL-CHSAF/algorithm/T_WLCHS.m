function [Data] =T_WLCHS(Input,noise,unknown_w,L,mu,lamada)

N = length(Input);
EMSE = zeros(1,N);
MSD = zeros(1,N);
[Rz,Pz,sigma_v,sigma_cv] = Sta_signal(Input,noise,L);
Q = unknown_w*unknown_w';
G = unknown_w*transpose(unknown_w);

for i = L:N
    
  A1  =2*Rz*Q*Rz*Q+Pz*conj(Q)*conj(Pz)*Q;
  A2  =Rz*Q*Rz + Pz*conj(Q)*conj(Pz) + trace(Rz*Q)*Rz;
  B1  =2*Rz*Q*Rz*G+Pz*conj(Q)*conj(Pz)*G;
  B2  =2*Rz*G*conj(Rz) + trace(conj(Pz)*G)*Pz;
  Q = Q -  mu*lamada*(Rz*Q+Q*Rz)+(5/6)*mu*lamada^3*(A1+conj(sigma_cv)*Pz*conj(G)+2*sigma_v*Rz*Q+A1'+sigma_cv*G*conj(Pz)+2*sigma_v*Q*Rz)...
      + mu^2*lamada^2*(A2+sigma_v*Rz);
  
  G = G -  mu*lamada*(Rz*G+G*conj(Rz))+(5/6)*mu*lamada^3*(B1+conj(sigma_cv)*Pz*conj(Q)+2*sigma_v*Rz*G+transpose(B1)+conj(sigma_cv)*Q*Pz+2*sigma_v*G*conj(Rz))...
      + mu^2*lamada^2*(B2+conj(sigma_cv)*Pz);
  
  MSD(i) = real(trace(Q));  
  EMSE(i) = real(trace(Q*Rz)); 

end
Data.TMSD  = MSD;
Data.TEMSE  = EMSE;
end





