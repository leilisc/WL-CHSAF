clear; 
addpath('./algorithm');
addpath('./signal'); 
% ¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
K = 100; 
N = 4000; 
MSDg = zeros(1,N); 
MSDb = zeros(1,N); 
MSDr = zeros(1,N); 
MSDk = zeros(1,N); 
MSDm = zeros(1,N); 
mu =0.1;
lamada = 0.5; 
alpha = 5;
%¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
for k = 1:K  
% ¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­     
 L = 8; 
unknown_w = unifrnd(-0.5,0.5,2*L,1) + 1j*unifrnd(-0.5,0.5,2*L,1);
unknown_w2 = -unknown_w;
% ¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­     
      u2 = sqrt(0.5)*(randn(1,N) + 1j*randn(1,N)); 
%      u2 = Noncircular_CG(N);

      Input = u2;   
% ¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­     
     D1 = filter(conj(unknown_w(1:L)), 1, Input(1:N/2)) + filter(conj(unknown_w(L+1:2*L)), 1, conj(Input(1:N/2))); %output
     D2 = filter(conj(unknown_w2(1:L)), 1, Input(N/2+1:N)) + filter(conj(unknown_w2(L+1:2*L)), 1, conj(Input(N/2+1:N))); %output
     D =[D1 D2];
     noise = Complex_CG(N,0.008,0.002,0.05,8,2);
     Output = D + noise;
     Lw =  length(unknown_w); 
%¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
   
    [Da]= WLCHS2(Input,Output,D,unknown_w,Lw,mu,0.2);
    MSDg = MSDg + Da.MSD;  
    [Da]= WLCHS2(Input,Output,D,unknown_w,Lw,mu,0.4);
    MSDb = MSDb + Da.MSD;
    [Da]= WLCHS2(Input,Output,D,unknown_w,Lw,mu,0.6);
    MSDr = MSDr + Da.MSD;
    [Da]= WLCHS2(Input,Output,D,unknown_w,Lw,mu,0.8);
    MSDk = MSDk + Da.MSD;
    [Da]= WLCHS2(Input,Output,D,unknown_w,Lw,mu,1);
    MSDm = MSDm + Da.MSD;
    
end

%¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
%¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
MSDg = MSDg/K; 
MSDb = MSDb/K;  
MSDr = MSDr/K;  
MSDk = MSDk/K;  
MSDm = MSDm/K;  

n = 1:N; 
figure(1);
plot(n,10*log10(MSDg(n)),'g','linewidth',1);
hold on;
plot(n,10*log10(MSDb(n)),'b','linewidth',1);
hold on;
plot(n,10*log10(MSDr(n)),'r','linewidth',1);
hold on;
plot(n,10*log10(MSDk(n)),'k','linewidth',1);
hold on;
plot(n,10*log10(MSDm(n)),'m','linewidth',1);
hold on;
xlabel('iteration','FontSize',15,'FontName','Times New Roman');
ylabel('MSD(dB)','FontSize',15,'FontName','Times New Roman');
h=legend('WL-CHSAF \lambda=0.2','WL-CHSAF \lambda=0.4','WL-CHSAF \lambda=0.6','WL-CHSAF \lambda=0.8','WL-CHSAF \lambda=1');
set(h,'FontSize',10,'FontName','Times New Roman');