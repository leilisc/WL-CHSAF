function [d] =ARMA_signal(N)
%Generate complex ARMA signal
noise = sqrt(0.5)*(randn(1,N) + 1j*randn(1,N));
d = zeros(1,N);
d(1) = noise(1);
 for k =2:N 
  d(k) = 0.6*d(k-1) + 2*noise(k) + 0.5*conj(noise(k)) + 0.1*noise(k-1) + 2*conj(noise(k-1));
 end
 
end

