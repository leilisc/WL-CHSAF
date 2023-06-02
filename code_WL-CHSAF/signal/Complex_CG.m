function Y = Complex_CG(N,sigma11,sigma12,P,sigma21,sigma22)
%Generate mixed complex noise
noise1 = sqrt(sigma11)*rand(1,N) + 1j*sqrt(sigma12)*rand(1,N);
noise2 =  binornd(1,P,1,N).*(sqrt(sigma21)*rand(1,N)) + 1j*binornd(1,P,1,N).*(sqrt(sigma22)*rand(1,N));
Y = noise1 + noise2;
end



