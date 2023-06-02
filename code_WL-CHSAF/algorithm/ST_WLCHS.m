function [Data] =ST_WLCHS(Input,noise,L,mu,lamada)

[Rz,Pz,sigma_v,sigma_cv] = Sta_signal(Input,noise,L);

Lw=2*L;
V_I = reshape(eye(Lw),Lw^2,1);
V_Rz  = reshape(Rz,Lw^2,1);
V_Rzt = reshape(transpose(Rz),Lw^2,1);
V_Pz  = reshape(Pz,Lw^2,1);
V_Pzc = reshape(conj(Pz),Lw^2,1);

 
 F1= eye(Lw^2) - mu*lamada*kron(eye(Lw),Rz) - mu*lamada*kron(conj(Rz),eye(Lw)) + 2*mu*sigma_v*lamada^3*(5/6)*(kron(eye(Lw),Rz)+kron(conj(Rz),eye(Lw))) + mu^2*lamada^2*(kron(conj(Rz),Rz)+V_Rz*transpose(V_Rzt)); 
 F2= mu^2*lamada^2*kron(conj(Pz),Pz);
 F3= mu*sigma_cv*lamada^3*(5/6)*kron(conj(Pz),eye(Lw));
 F4= mu*conj(sigma_cv)*lamada^3*(5/6)*kron(eye(Lw),Pz);
 H1= eye(Lw^2) - mu*lamada*kron(eye(Lw),Rz) - mu*lamada*kron(Rz,eye(Lw)) + 2*mu*sigma_v*lamada^3*(5/6)*(kron(eye(Lw),Rz)+kron(Rz,eye(Lw))) + mu^2*lamada^2*(2*kron(Rz,Rz)+V_Pz*transpose(V_Pzc)); 
 
 J = eye(Lw^2)-F1-F3*(eye(Lw^2)-H1)^(-1)*conj(F3)-F4*(eye(Lw^2)-conj(H1))^(-1)*conj(F4);
 T = F2+F3*(eye(Lw^2)-H1)^(-1)*F4 + F4*(eye(Lw^2)-conj(H1))^(-1)*F3;
 M = F3*(eye(Lw^2)-H1)^(-1)*mu^2*lamada^2*conj(sigma_cv)*V_Pz + F4*(eye(Lw^2)-conj(H1))^(-1)*mu^2*lamada^2*sigma_cv*V_Pzc + mu^2*lamada^2*sigma_v*V_Rz;
 
 V_Q = (J-T*conj(J^(-1))*conj(T))^(-1)*(T*conj(J^(-1))*conj(M)+M);
 
 Data.STMSD = real(transpose(V_I)*V_Q);  
 Data.STEMSE = real(transpose(V_Rzt)*V_Q); 
end

