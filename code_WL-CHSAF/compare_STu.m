clc; 
clear; 
addpath('./algorithm'); 
addpath('./signal'); 
% ¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
K = 50; 
N = 5000; 

lamada = 0.2;
mu=(0.05:0.02:0.4);
%¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
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

 L = 4; 
 unknown_w = unifrnd(-0.5,0.5,2*L,1) + 1j*unifrnd(-0.5,0.5,2*L,1);

% ¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­     
       Input1 = sqrt(0.5)*randn(1,N) + sqrt(0.1)*1j*randn(1,N); 
       Input2 = sqrt(0.7)*randn(1,N) + sqrt(0.3)*1j*randn(1,N); 
       Input3 = sqrt(0.9)*randn(1,N) + sqrt(0.5)*1j*randn(1,N); 

%        Input1 = sqrt(0.5)*randn(1,N) + sqrt(0.5)*1j*randn(1,N); 
%        Input2 = sqrt(0.7)*randn(1,N) + sqrt(0.3)*1j*randn(1,N); 
%        Input3 = sqrt(0.9)*randn(1,N) + sqrt(0.1)*1j*randn(1,N); 

% ¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­     
     D1 = filter(conj(unknown_w(1:L)), 1, Input1) + filter(conj(unknown_w(L+1:2*L)), 1, conj(Input1)); %output
     D2 = filter(conj(unknown_w(1:L)), 1, Input2) + filter(conj(unknown_w(L+1:2*L)), 1, conj(Input2)); %output
     D3 = filter(conj(unknown_w(1:L)), 1, Input3) + filter(conj(unknown_w(L+1:2*L)), 1, conj(Input3)); %output
     noise = Complex_CG(N,0.008,0.002,0,90,10);
     Output1 = D1 + noise;
     Output2 = D2 + noise;
     Output3 = D3 + noise;
     Lw =  length(unknown_w); 
%¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
   
    [Da]= WLCHS(Input1,Output1,D1,unknown_w,Lw,mu(n),lamada);
    MSDg = MSDg + Da.MSD;  
    [Da]= WLCHS(Input2,Output2,D2,unknown_w,Lw,mu(n),lamada);
    MSDb = MSDb + Da.MSD;
    [Da]= WLCHS(Input3,Output3,D3,unknown_w,Lw,mu(n),lamada);
    MSDr = MSDr + Da.MSD;


    [Da]= ST_WLCHS(Input1,noise,L,mu(n),lamada);
    TMSDg = TMSDg + Da.STMSD;  
    [Da]= ST_WLCHS(Input2,noise,L,mu(n),lamada);
    TMSDb = TMSDb + Da.STMSD;  
    [Da]= ST_WLCHS(Input3,noise,L,mu(n),lamada);
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
%¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
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
h=legend('Simulation ($\sigma_{u}^2=0.6$ $\tilde{\sigma}_{u}^{2}=0.4$)','Simulation ($\sigma_{u}^2=1.0$  $\tilde{\sigma}_{u}^{2}=0.4$)','Simulation ($\sigma_{u}^2=1.4$ $\tilde{\sigma}_{u}^{2}=0.4$)','Theory ($\sigma_{u}^2=0.6$ $\tilde{\sigma}_{u}^{2}=0.4$)','Theory ($\sigma_{u}^2=1.0$ $\tilde{\sigma}_{u}^{2}=0.4$)','Theory ($\sigma_{u}^2=1.4$ $\tilde{\sigma}_{u}^{2}=0.4$)');
% h=legend('Simulation ($\sigma_{u}^2=1$ $\tilde{\sigma}_{u}^{2}=0$)','Simulation ($\sigma_{u}^2=1$  $\tilde{\sigma}_{u}^{2}=0.4$)','Simulation ($\sigma_{u}^2=1$ $\tilde{\sigma}_{u}^{2}=0.8$)','Theory ($\sigma_{u}^2=1$ $\tilde{\sigma}_{u}^{2}=0$)','Theory ($\sigma_{u}^2=1$ $\tilde{\sigma}_{u}^{2}=0.4$)','Theory ($\sigma_{u}^2=1$ $\tilde{\sigma}_{u}^{2}=0.8$)');
set(h,'FontSize',10,'FontName','Times New Roman', 'interpreter', 'latex');