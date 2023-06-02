clear; 
e = 0:0.001:10; 
lamada = 0.1;

Y1 = g_F(e,0.2);
Y2 = g_F(e,0.4);
Y3 = g_F(e,0.6);
Y4 = g_F(e,0.8);
Y5 = g_F(e,1);

figure(1);
plot(e,Y1,'g','linewidth',1.5);
hold on;
plot(e,Y2,'b','linewidth',1.5);
hold on;
plot(e,Y3,'r','linewidth',1.5);
hold on;
plot(e,Y4,'m','linewidth',1.5);
hold on;
plot(e,Y5,'y','linewidth',1.5);
hold on;
xlabel('|e(k)|','FontSize',15,'FontName','Times New Roman');
% ylabel('sech(\lambda|e(k)|)','FontSize',15,'FontName','Times New Roman');
ylabel('q(|e(k)|)','FontSize',15,'FontName','Times New Roman');
h=legend('\lambda=0.2','\lambda=0.4','\lambda=0.6','\lambda=0.8','\lambda=1');
set(h,'FontSize',10,'FontName','Times New Roman');



 function [Y]= g_F(e,lamada)
    L=length(e);
    Y=zeros(1,L);
    for i =1:L
%       Y(i)=sech(lamada*abs(e(i)));  
      Y(i)=(sech(lamada*abs(e(i))))*tanh(lamada*abs(e(i)))/abs(e(i));
    end
 end