clc; 
clear; 
addpath('./algorithm'); 
addpath('./signal'); 
% ¡­¡­¡­¡­¡­¡­¡­¡­parameter¡­¡­¡­¡­¡­¡­¡­¡­¡­
K = 500; 
N = 2000; 
MSDg = zeros(1,N); 
MSDb = zeros(1,N); 
MSDr = zeros(1,N); 
MSDk = zeros(1,N); 
MSDm = zeros(1,N); 
MSDy = zeros(1,N); 

lamada = 0.5;
%¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
for k = 1:K  
% ¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
L = 4; 
unknown_w = unifrnd(-0.5,0.5,L,1) + 1j*unifrnd(-0.5,0.5,L,1);

% ¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­     
   Input = AR1_signal(N);
%              Input = ARMA_signal(N);
   
  D = [Input(2:end) 0];
  Output = [Input(2:end) 0];
unknown_w2= [unknown_w;zeros(L,1)];
Lw = length(unknown_w2); % 
%¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
   
    [Da]= CHS(Input,Output,D,unknown_w2,Lw,0.005,lamada);
    MSDg = MSDg + Da.MSE;  
    [Da]= WLCHS(Input,Output,D,unknown_w2,Lw,0.005,lamada);
    MSDb = MSDb + Da.MSE;
    [Da]= CHS(Input,Output,D,unknown_w2,Lw,0.02,lamada);
    MSDr = MSDr + Da.MSE;
    [Da]= WLCHS(Input,Output,D,unknown_w2,Lw,0.02,lamada);
    MSDk = MSDk + Da.MSE;
 
end

%¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
%¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
MSDg = MSDg/K;  
MSDb = MSDb/K;  
MSDr = MSDr/K;  
MSDk = MSDk/K;  
MSDm = MSDm/K;  
MSDy = MSDy/K;  



figure(1);
plot(10*log10(smooth(MSDg,10)),'g','linewidth',1);
hold on;
plot(10*log10(smooth(MSDb,10)),'b','linewidth',1);
hold on;
plot(10*log10(smooth(MSDr,10)),'r','linewidth',1);
hold on;
plot(10*log10(smooth(MSDk,10)),'k','linewidth',1);
hold on;
xlabel('iteration','FontSize',15,'FontName','Times New Roman');
ylabel('MSE(dB)','FontSize',15,'FontName','Times New Roman');
h=legend('CHSAF \mu=0.005','WL-CHSAF \mu=0.005 ','CHSAF \mu=0.02','WL-CHSAF \mu=0.02');
% h=legend('CHSAF \mu=0.005','WL-CHSAF \mu=0.005 ','CHSAF \mu=0.01','WL-CHSAF \mu=0.01');
set(h,'FontSize',10,'FontName','Times New Roman');