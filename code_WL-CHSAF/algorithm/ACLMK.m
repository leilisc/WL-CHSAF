function [Data] = ACLMK(Input,Output,D,unknown_w,Lw,mu)
% ACLMK algorithm
  L= Lw/2;
  N = length(Input);
  Y = zeros(1,N); 
  e = zeros(1,N);
  ea = zeros(1,N); 
  MSD = zeros(1,N); 
  h = zeros(L,1); 
  g = zeros(L,1);
  W = [h;g]; 
  lambda = 0.995;
  P = 0;
  for i = L:N
      X = transpose(Input(i:-1:i-L+1)); 
      Y(i) = h'*X + g'*conj(X);
      e(i) = Output(i) - Y(i); 
      ea(i) = D(i) - Y(i);
      P = lambda*P + e(i)*conj(e(i));
      k1=(3*P-e(i)*conj(e(i)))*conj(e(i));
      h = h + mu*k1*X;
      g = g + mu*k1*conj(X);
      W = [h;g];
      MSD(i) = (norm(W - unknown_w))^2;
  end    
Data.MSD  = MSD;
Data.MSE  = abs(e).^2;
Data.EMSE = abs(ea).^2;
end