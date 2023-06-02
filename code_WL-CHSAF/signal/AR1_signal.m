function [d] = AR1_signal(N)
%Generate complex AR1 signal
noise = sqrt(0.5)*randn(1,N) + sqrt(0.5)*1j*randn(1,N);
d = zeros(1,N);
d(1) = noise(1);
 for k = 2:N 
  d(k) = 0.8*d(k-1) + noise(k);
 end
end

