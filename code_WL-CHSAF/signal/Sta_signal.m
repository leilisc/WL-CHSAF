function [Rz,Pz,sigma_v,sigma_cv] = Sta_signal(Input,noise,L)

N = length(Input);
Lw=2*L;
Rz = zeros(Lw,Lw);
Pz = zeros(Lw,Lw);
sigma_v = 0;
sigma_cv= 0;
 for i = L:N
      X = transpose(Input(i:-1:i-L+1)); 
      Z = [X; conj(X)]; 
      Rz =  Rz+Z*Z';   
      Pz =  Pz+Z*transpose(Z);   
      sigma_v= sigma_v + noise(i)*conj(noise(i));
      sigma_cv= sigma_cv + noise(i)^2;
 end
     Rz =  Rz/(N-L+1);
     Pz =  Pz/(N-L+1);
     sigma_v =  sigma_v/(N-L+1);
     sigma_cv =  sigma_cv/(N-L+1);
end

