%%
%WCMC-Q -Research Dept.
% Date: 21/03/2011
% Author: Susanne Awad
function outy=Risk_HIV_Model_Changed(param)
global mu N0 epsilon nst rst nalpha h fb NST date p n tau B
%define the number of stages and risk groups.
nst =4;% str2num(answer2{1});
rst = 10;%str2num(answer{1});
prec=0.2;%.2;
%%
%Time scale
t0=1970;         %Start time
tf=2030;         %Stop time
dt=0.5;         %time interval
tspan=t0:dt:tf;  %timespan to use in the ode45
d=(date-t0)*2+1;
d1=(1990-t0)*2+1;
d2=(2009-t0)*2+1;
%% Biological co-factor
NST=2*nst;
fb=0;%0.5;
h=1;%0.7;
%%
% The Demographic Parameters

N0tot=1000;           % the total initial population
% N0=f.*N0tot;         % the fraction of the inital population in THE N risk group
m=1:rst;
f = gampdf(m,1.9,2);%exppdf(m,4.4325);%exp(-0.001*m);% m=1:.935:10;x=1;for b=1:10000
% f(:,b) = gampdf(m,x,2.4);
% x=x+1e-3;
% end
% an=sum(f,1);
% plot(an)
% 1+1000*1e-3
% 
% ans =
% 
%      2

N0=(f.*N0tot);      
mu=1/35;             %The duration of the active sexual lifetime (15-49yr)
%%
% Transmission Parameters
 p(2)=0.03604;         %unifrnd(0.03604*(1-prec),0.03604*(1+prec));pb of transmission-acute stage in r risk groups
 p(3)=0.0008;         %unifrnd(0.0008*(1-prec),0.0008*(1+prec));%pb of transmission-latent stage in r risk groups
 p(4)=0.0042;         %unifrnd(0.0042*(1-prec),0.0042*(1+prec));%pb of transmission-advanced stage in r risk groups
p(6)=p(2);
p(7)=p(3);
p(8)=p(4);

%Progression parameters
 omg(2)=365/49;       %unifrnd(365/49*(1-prec),365/49*(1+prec));%rate at which people move from acute to latent stage,2 month
 omg(3)=1/7.59;          %unifrnd(1/7.59*(1-prec),1/7.59*(1+prec));%rate at which people move from the latent to advanced stage, 10 yrs
 omg(4)=1/2;%unifrnd(1/2*(1-prec),1/2*(1+prec));%rate at which people move from the latent to advanced stage, 2 yrs
omg(6)=omg(2);
omg(7)=omg(3);
omg(8)=omg(4);

% Sexual behavior parameters
 n(2)=10.6;           %unifrnd(10.6*(1-prec),10.6*(1+prec));%coital frequency per month-acute stage
 n(3)=11;             %unifrnd(11*(1-prec),11*(1+prec));%coital frequency per month-latent stage
 n(4)=7.1;            %unifrnd(7.1*(1-prec),7.1*(1+prec));%coital frequency per month-advanced stage
n(6)=n(2);
n(7)=n(3);
n(8)=n(4);

tau(1:rst,1:rst)=6;  %unifrnd(6*(1-prec),6*(1+prec));%duration of partnership=6months

%%
%Epsilon for the G-matrix
epsilon =0.3;%unifrnd(0.3*(1-prec),0.3*(1+prec));% str2num(answer{1});
%%
%The transmission probability per parnetship
q1=zeros(nst,rst,rst);
for alpha=2:nst        %3 stages of HIV
for j = 1:rst  
    for i= 1:rst
        q1(alpha,j,i) =1-(1-p(alpha))^(n(alpha)*tau(j,i));
    end
end
end
%FOR Intervention
q2=zeros(nst,rst,rst);
for alpha=nst+2:NST        %3 stages of HIV
for j = 1:rst  
    for i= 1:rst
        q2(alpha,j,i) =1-(1-p(alpha))^(n(alpha)*tau(j,i));
    end
end
end
%%
%The new sexual partner acquisition rate
Z=param(2);
zeta_s=param(3);
zeta_y=param(4);
A0=param(1);
 B=1.7;%unifrnd(1.7*(1-prec),1.7*(1+prec));%
for i=1:rst
rho_unif(i)=A0*i^B;
end

%%
%Initial conditions
 initprev=param(5);%0.01/100; %initial prevalence for all risk groups at the initial conditions
x0L=zeros(rst,NST);
x0L(:,1)= (1-fb)*N0(:)-initprev*(1-fb)*N0(:);
x0L(:,2)= initprev*(1-fb)*N0(:);
x0L(:,5)=(fb)*N0(:)-initprev*(fb)*N0(:);
x0L(:,6)=initprev*(fb)*N0(:);
x0=zeros(NST*rst,1);
i=1:nst;
for j=1:rst
    isj= i+(j-1)*NST;
    x0(isj)=x0L(j,i);
end
%%
%Solve the differential equation for the HIV system of equations model
[t,x]=ode45(@risk_HIV,tspan,x0,[],omg,q1,q2,rho_unif,zeta_s,zeta_y,Z);       
%%
% [dxL,lmp]=risk_HIV(t,x,omg,rho,q1,q2);
 runs=length(t);

% %Total Population
Ntot=sum(x,2);
% 
% %%
Y = zeros(runs,rst);
i=2;

for j=1:rst
   for tx=1:runs
       Y(tx,j)=sum(x(tx,i:i+(nst-2)));
 
end
i=i+NST;
end
Ytot=sum(Y,2);
% %%

%Total Prevalence
PrevTOT= Ytot./Ntot*100;
% figure(1)
plot(t,PrevTOT)
outy=PrevTOT(d1:2:d2,1)';
end
