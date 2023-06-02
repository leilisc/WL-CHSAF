function Y = Noncircular_CG(N)
%Generate non-circular complex Gaussian signal£ºvariance=1,complementary variance=0.4+0.3*j.
Xr = sqrt(0.7)*randn(1,N) ;
Xi = (3/14)*Xr + sqrt(105/392)*randn(1,N);
Y = Xr + 1j*Xi;
end



