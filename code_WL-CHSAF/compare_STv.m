clc; 
clear; 
addpath('./algorithm'); 
addpath('./signal'); 
% ¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
K = 50; 
N = 5000; 

lamada = 0.2;
mu=(0.05:0.02:0.4);
%¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
for n = 1:length(mu) 
    
MSDg = zeros(1,N); 
MSDb = zeros(1,N); 
MSDr = zeros(1,N); 
MSDk = zeros(1,N); 

TMSDg = 0; 
TMSDb = 0; 
TMSDr = 0; 
TMSDk = 0; 
for k = 1:K  
% ¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­     
L = 4; 
unknown_w = unifrnd(-0.5,0.5,2*L,1) + 1j*unifrnd(-0.5,0.5,2*L,1);

    Input = sqrt(0.5)*(randn(1,N) + 1j*randn(1,N));
% ¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­     
     D = filter(conj(unknown_w(1:L)), 1, Input) + filter(conj(unknown_w(L+1:2*L)), 1, conj(Input)); 
     noise1 = Complex_CG(N,0.005,0.001,0,90,10);
     noise2 = Complex_CG(N,0.007,0.003,0,90,10);
     noise3 = Complex_CG(N,0.009,0.005,0,90,10);
     
%      noise1 = Complex_CG(N,0.01,0.01,0,90,10);
%      noise2 = Complex_CG(N,0.012,0.008,0,90,10);
%      noise3 = Complex_CG(N,0.014,0.006,0,90,10);
     
     Output1 = D + noise1;
     Output2 = D + noise2;
     Output3 = D + noise3;
     Lw =  length(unknown_w); 
%¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
   
    [Da]= WLCHS(Input,Output1,D,unknown_w,Lw,mu(n),lamada);
    MSDg = MSDg + Da.MSD;  
    [Da]= WLCHS(Input,Output2,D,unknown_w,Lw,mu(n),lamada);
    MSDb = MSDb + Da.MSD;
    [Da]= WLCHS(Input,Output3,D,unknown_w,Lw,mu(n),lamada);
    MSDr = MSDr + Da.MSD;


    [Da]= ST_WLCHS(Input,noise1,L,mu(n),lamada);
    TMSDg = TMSDg + Da.STMSD;  
    [Da]= ST_WLCHS(Input,noise2,L,mu(n),lamada);
    TMSDb = TMSDb + Da.STMSD;  
    [Da]= ST_WLCHS(Input,noise3,L,mu(n),lamada);
    TMSDr = TMSDr + Da.STMSD;  

       
end


SMSDg(n) = mean(MSDg(end-500:end))/K;
SMSDb(n) = mean(MSDb(end-500:end))/K;
SMSDr(n) = mean(MSDr(end-500:end))/K;
SMSDk(n) = mean(MSDk(end-500:end))/K;


STMSDg(n) = TMSDg/K;
STMSDb(n) = TMSDb/K;
STMSDr(n) = TMSDr/K;
STMSDk(n) = TMSDk/K;

end
%¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­

figure(1);
plot(mu,10*log10(SMSDg),'--go','linewidth',1);
hold on;
plot(mu,10*log10(SMSDb),'--bo','linewidth',1);
hold on;
plot(mu,10*log10(SMSDr),'--ro','linewidth',1);
hold on;

plot(mu,10*log10(STMSDg),'-g^','linewidth',1);
hold on;
plot(mu,10*log10(STMSDb),'-b^','linewidth',1);
hold on;
plot(mu,10*log10(STMSDr),'-r^','linewidth',1);
hold on;

xlabel('\mu ','FontSize',15,'FontName','Times New Roman');
ylabel('Steady-state MSD(dB)','FontSize',15,'FontName','Times New Roman');
h=legend('Simulation ($\sigma_{v}^2=0.006$ $\tilde{\sigma}_{v}^{2}=0.004$)','Simulation ($\sigma_{v}^2=0.01$  $\tilde{\sigma}_{v}^{2}=0.004$)','Simulation ($\sigma_{v}^2=0.014$ $\tilde{\sigma}_{v}^{2}=0.004$)','Theory ($\sigma_{v}^2=0.006$ $\tilde{\sigma}_{v}^{2}=0.004$)','Theory ($\sigma_{v}^2=0.01$ $\tilde{\sigma}_{v}^{2}=0.004$)','Theory ($\sigma_{v}^2=0.014$ $\tilde{\sigma}_{v}^{2}=0.004$)');
% h=legend('Simulation ($\sigma_{v}^2=0.02$ $\tilde{\sigma}_{v}^{2}=0$)','Simulation ($\sigma_{v}^2=0.02$  $\tilde{\sigma}_{v}^{2}=0.004$)','Simulation ($\sigma_{v}^2=0.02$ $\tilde{\sigma}_{v}^{2}=0.008$)','Theory ($\sigma_{v}^2=0.02$ $\tilde{\sigma}_{v}^{2}=0$)','Theory ($\sigma_{v}^2=0.02$ $\tilde{\sigma}_{u}^{2}=0.004$)','Theory ($\sigma_{v}^2=0.02$ $\tilde{\sigma}_{v}^{2}=0.008$)');
set(h,'FontSize',10,'FontName','Times New Roman', 'interpreter', 'latex');