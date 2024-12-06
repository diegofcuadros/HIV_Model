function [dxL,lambdaTot]=risk_HIV(t,x,omg,q1,q2,rho_unif,zeta_s,zeta_y,Z)
% Right hand sides of system of D.E.'s
%  dx(1)/dt=...
%  dx(2)/dt=...
%  dx(3)/dt=...
%  dx(4)/dt=... where S=x(1),I_1=x(2), I_2=x(3), I_3=x(4).
% Both x and xprime are column vectors.
global mu N0 epsilon  nst rst h fb NST p n tau B

%xprime=zeros(nst*rst,1);
xt=zeros(rst,NST);
i=1:NST;
for j=1:rst
    isj= i+(j-1)*NST;
    xt(j,i)=x(isj);
end
Risk=0;
for i=1:NST*rst
    Risk(i)= x(i);
end
  rho=rhot(t,rho_unif,zeta_s,zeta_y,Z);%zeros(NST,rst);
% hold;plot(t,rho(1,4));hold
%%
%recall the force of infection for the risk groups
lambdaTot=lambda(q1,q2,rho,Risk);
%%
%the Differential Equations
dx=zeros(rst,NST);
for i=1:rst
dx(i,1)= (1-fb)*mu*N0(i)- mu*xt(i,1)-lambdaTot(i)*xt(i,1);
dx(i,2)= lambdaTot(i)*xt(i,1) - mu*xt(i,2)-omg(2)*xt(i,2);
for j=3:nst
dx(i,j)= omg(j-1)*xt(i,j-1)-mu*xt(i,j)-omg(j)*xt(i,j);
end
dx(i,5)= fb*mu*N0(i)- mu*xt(i,5)-h*lambdaTot(i)*xt(i,5);
dx(i,6)= h*lambdaTot(i)*xt(i,5) - mu*xt(i,6)-omg(6)*xt(i,6);
dx(i,7)= omg(6)*xt(i,6)-mu*xt(i,7)-omg(7)*xt(i,7);
dx(i,8)= omg(7)*xt(i,7)-mu*xt(i,8)-omg(8)*xt(i,8);
end
dxL=zeros(rst*NST,1);
i=1:NST;
for j=1:rst
    isj=i+(j-1)*NST;
dxL(isj)= dx(j,i);
end
end

