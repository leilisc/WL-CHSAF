function [Data] = CHS(Input,Output,D,unknown_w,L,mu,lamada)
% CHS algorithm
  N = length(Input);
  Y = zeros(1,N); 
  e = zeros(1,N); 
  ea = zeros(1,N); 
  MSD = zeros(1,N); 
  W = zeros(L,1); 
  for i = L:N
     X = transpose(Input(i:-1:i-L+1)); 
     Y(i) = W'*X;
      e(i) = Output(i) - Y(i); 
      ea(i) = D(i) - Y(i); 
      fai=sech(lamada*abs(e(i)))*tanh(lamada*abs(e(i)))*sign(conj(e(i)));
      W = W + mu*fai*X; 
      MSD(i) = (norm(W - unknown_w))^2;
  end    
Data.MSD  = MSD;
Data.MSE  = abs(e).^2;
Data.EMSE = abs(ea).^2;
end