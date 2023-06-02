function [Data] =ACILSE(Input,Output,D,unknown_w,Lw,mu,beta)
% ACILSE algorithm
  L= Lw/2;
  N = length(Input);
  Y = zeros(1,N); 
  e = zeros(1,N);
  ea = zeros(1,N); 
  lambda = 0.995;
  delta1 = 0;
  delta2 = 0;
  MSD = zeros(1,N); 
  h = zeros(L,1); 
  g = zeros(L,1);
  W = [h;g];  
  for i = L:N
     X = transpose(Input(i:-1:i-L+1)); 
     Y(i) = h'*X + g'*conj(X);
     e(i) = Output(i) - Y(i); 
     ea(i) = D(i) - Y(i); 
     delta1 = lambda*delta1 + (1-lambda)*e(i)*conj(e(i));
     delta2 = lambda*delta2 + (1-lambda)*e(i)*e(i);
     pho=delta2/delta1;
     theta=exp(-beta*delta1);
     fai = conj(e(i))-theta*conj(pho)*e(i);
     h = h + mu*fai*X;
     g = g + mu*fai*conj(X);
     W = [h;g];
     MSD(i) = (norm(W - unknown_w))^2;
  end
Data.MSD  = MSD;
Data.MSE  = abs(e).^2;
Data.EMSE = abs(ea).^2;
end