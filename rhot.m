function rho=rhot(t,rho_unif,zeta_s,zeta_y,Z)
global mu N0 epsilon  nst rst h fb NST
rho(1:NST,1:rst)=0;
for i=1:rst
    rho(1:NST,i)=rho_unif(i)+(Z*rho_unif(i))./(1+exp((t-zeta_s)./zeta_y));
end

end