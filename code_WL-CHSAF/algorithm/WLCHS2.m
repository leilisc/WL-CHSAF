function [Data] = WLCHS2(Input,Output,D,unknown_w,Lw,mu,lamada)
% WLCHS algorithm
  L= Lw/2;
  N = length(Input);
  Y = zeros(1,N); 
  e = zeros(1,N);
  ea = zeros(1,N); 
  MSD = zeros(1,N); 
  h = zeros(L,1); 
  g = zeros(L,1);
  W = [h;g]; 
  for i = L:N
     X = transpose(Input(i:-1:i-L+1)); 
     Y(i) = h'*X + g'*conj(X);
     e(i) = Output(i) - Y(i); 
     ea(i) = D(i) - Y(i); 
     fai=sech(lamada*abs(e(i)))*tanh(lamada*abs(e(i)))*sign(conj(e(i)));
     h = h + mu*fai*X;
     g = g + mu*fai*conj(X);
     W = [h;g];
      if i <= N/2
           MSD(i) =  (norm(W - unknown_w))^2;
         else
           MSD(i) =  (norm(W + unknown_w))^2;
      end
  end
Data.MSD  = MSD;
Data.MSE  = abs(e).^2;
Data.EMSE = abs(ea).^2;
end