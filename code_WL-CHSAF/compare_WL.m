
clear; 
addpath('./algorithm');
addpath('./signal'); 
% ¡­¡­¡­¡­¡­¡­¡­¡­parameter¡­¡­¡­¡­¡­¡­¡­¡­¡­
K = 100; 
N = 3000; 
MSDg = zeros(1,N); 
MSDb = zeros(1,N); 
MSDr = zeros(1,N); 
MSDk = zeros(1,N); 
MSDm = zeros(1,N); 
MSDy = zeros(1,N); 

lamada = 0.5;
alpha = 5;
%¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
for k = 1:K  
 
L = 8; 
unknown_w = unifrnd(-0.5,0.5,2*L,1) + 1j*unifrnd(-0.5,0.5,2*L,1);
% ¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­     
      u2 = sqrt(0.5)*(randn(1,N) + 1j*randn(1,N)); 
%      u2 = Noncircular_CG(N);
      Input = u2;   
% ¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­     
     D = filter(conj(unknown_w(1:L)), 1, Input) + filter(conj(unknown_w(L+1:2*L)), 1, conj(Input)); %output
     noise = Complex_CG(N,0.008,0.002,0.01,8,2);
     Output = D + noise;
     Lw =  length(unknown_w); 
 
%¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
   

    [Da]= CHS(Input,Output,D,unknown_w,Lw,0.02,lamada);
    MSDg = MSDg + Da.MSD;
    [Da]= WLCHS(Input,Output,D,unknown_w,Lw,0.02,lamada);
    MSDb = MSDb + Da.MSD;
    [Da]= CHS(Input,Output,D,unknown_w,Lw,0.05,lamada);
    MSDr = MSDr + Da.MSD;
    [Da]= WLCHS(Input,Output,D,unknown_w,Lw,0.05,lamada);
    MSDk = MSDk + Da.MSD;
    [Da]= CHS(Input,Output,D,unknown_w,Lw,0.1,lamada);
    MSDm = MSDm + Da.MSD;
    [Da]= WLCHS(Input,Output,D,unknown_w,Lw,0.1,lamada);
    MSDy = MSDy + Da.MSD;

  
end

%¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
%¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
MSDg = MSDg/K;  
MSDb = MSDb/K;  
MSDr = MSDr/K;  
MSDk = MSDk/K;  
MSDm = MSDm/K;  
MSDy = MSDy/K;  

n = 1:1:N; 
figure(1);
plot(n,10*log10(MSDg(n)),'g--','linewidth',1);
hold on;
plot(n,10*log10(MSDb(n)),'g','linewidth',1);
hold on;
plot(n,10*log10(MSDr(n)),'b--','linewidth',1);
hold on;
plot(n,10*log10(MSDk(n)),'b','linewidth',1);
hold on;
plot(n,10*log10(MSDm(n)),'r--','linewidth',1);
hold on;
plot(n,10*log10(MSDy(n)),'r','linewidth',1);
hold on;

xlabel('iteration','FontSize',15,'FontName','Times New Roman');
ylabel('MSD(dB)','FontSize',15,'FontName','Times New Roman');
h=legend('CHSAF \mu=0.02','WL-CHSAF \mu=0.02','CHSAF \mu=0.05','WL-CHSAF \mu=0.05 ','CHASF \mu=0.1','WL-CHSAF \mu=0.1');
set(h,'FontSize',10,'FontName','Times New Roman');