function [Data] = MICCC(Input,Output,D,unknown_w,Lw,mu,sigma,pho)
% MICCC algorithm
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
     k1 = exp(-(abs(e(i))^2-real(pho*(conj(e(i)))^2))/(sigma^2*(1-abs(pho)^2)))/(pi*sigma^2*sqrt(1-abs(pho)^2));
     k2 =(conj(e(i))-conj(pho)*e(i))/(sigma^2*(1-abs(pho)^2));
     h = h + mu*k1*k2*X;
     g = g + mu*k1*k2*conj(X);
     W = [h;g];
     MSD(i) = (norm(W - unknown_w))^2;
  end
Data.MSD  = MSD;
Data.MSE  = abs(e).^2;
Data.EMSE = abs(ea).^2;
end