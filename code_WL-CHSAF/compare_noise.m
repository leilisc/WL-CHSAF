clear; 
addpath('./algorithm'); 
addpath('./signal'); 

% °≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠
K = 100; 
N = 3000; 
MSDg = zeros(1,N); 
MSDb = zeros(1,N); 
MSDr = zeros(1,N); 
MSDk = zeros(1,N); 
MSDm = zeros(1,N); 

%°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠
for k = 1:K  
% °≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠
L = 8; 
unknown_w = unifrnd(-0.5,0.5,2*L,1) + 1j*unifrnd(-0.5,0.5,2*L,1);
 
%       u2 = sqrt(0.5)*(randn(1,N) + 1j*randn(1,N)); 
     Input = Noncircular_CG(N);

     D = filter(conj(unknown_w(1:L)), 1, Input) + filter(conj(unknown_w(L+1:2*L)), 1, conj(Input)); %output
     noise = Complex_CG(N,0.008,0.002,0.01,8,2);
     Output = D + noise;
     Lw =  length(unknown_w); % filter oder
%°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠À„∑®°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠
   
    [Da]= ACLMS(Input,Output,D,unknown_w,Lw,0.02);
    MSDg = MSDg + Da.MSD;  
    [Da]= ACLMK(Input,Output,D,unknown_w,Lw,0.0001);
    MSDb = MSDb + Da.MSD;
    [Da]= ACILSE(Input,Output,D,unknown_w,Lw,0.02,2);
    MSDr = MSDr + Da.MSD;
    [Da]= MICCC(Input,Output,D,unknown_w,Lw,1.5,2,0.2);
    MSDk = MSDk + Da.MSD;
    [Da]= WLCHS(Input,Output,D,unknown_w,Lw,0.04,0.8);
    MSDm = MSDm + Da.MSD;
    
%     [Da]= ACLMS(Input,Output,D,unknown_w,Lw,0.02);
%     MSDg = MSDg + Da.MSD;  
%     [Da]= ACLMK(Input,Output,D,unknown_w,Lw,0.00004);
%     MSDb = MSDb + Da.MSD;
%     [Da]= ACILSE(Input,Output,D,unknown_w,Lw,0.02,2);
%     MSDr = MSDr + Da.MSD;
%     [Da]= MICCC(Input,Output,D,unknown_w,Lw,1.5,2,0.2);
%     MSDk = MSDk + Da.MSD;
%     [Da]= WLCHS(Input,Output,D,unknown_w,Lw,0.05,0.8);
%     MSDm = MSDm + Da.MSD;
end

%°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠output°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠°≠
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
h=legend('ACLMS \mu=0.02','ACLMK \mu=0.0001','WL-ILSE \mu=0.02 \beta=2','MICCC \mu=1.5 \sigma=2 \rho=0.2','WL-CHSAF \mu=0.04 \lambda=0.8');
% h=legend('ACLMS \mu=0.02','ACLMK \mu=0.00004','WL-ILSE \mu=0.02 \beta=2','MICCC \mu=1.5 \sigma=2 \rho=0.2','WL-CHSAF \mu=0.05 \lambda=0.8');
set(h,'FontSize',10,'FontName','Times New Roman');