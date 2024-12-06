rst=10;
des={ 'Sierra Leonne' 'Mali' 'Liberia' 'Ghana' 'Cote dIvoire' 'Cameroon' 'Tanzania' 'Kenya' 'Malawi' 'Zambia' 'Zimbabwe' 'Swaziland' 'Lesotho'};  %COUNTRIES AS A STRING
Data1=xlsread('Data','sheet1');
load('Final Estimates_all_new_1-8_new.mat')
  for i=1:length(Data1(1,:))
     data1=Data1(2:end-1,i)'; %D
  date=Data1(end,i);
 subplot(6,2,i)
 model_end=Risk_HIV_Model_Changed(Final_Estimatesb2(i,:));
 hold
  plot(1990:2009, Data1(2:21,i),'*r')
 title((des(1,i+1)));
 axis([1980 2020 0 (max(model_end)+5)])
 end
B=1.7;
%%
A0=Final_Estimatesb2(:,1);
% B=Final_Estimatesb2(:,2);
Z=Final_Estimatesb2(:,2);
zeta_s=Final_Estimatesb2(:,3);
zeta_y=Final_Estimatesb2(:,4);
t=1970:0.5:2030;
    for i=1:10
        for j=1:12
rho_unif(j,i)=A0(j,1).*i.^B;
    end
    end
for i=1:10
    for j=1:12
    for n=1:121 
rho(n,i,j)=rho_unif(j,i)+((Z(j,1)*rho_unif(j,i))./(1+exp((t(n)-zeta_s(j,1))./zeta_y(j,1))));
    end
    end
end
pd=(rho(1,10,:)-rho(end,10,:))./rho(1,10,:).*100;
pd(1,:)=pd(1,1,:);
figure(2)
subplot(3,1,1);plot(1:12,pd(1,:),'*r')
subplot(3,1,2);plot(1:12,zeta_y(:,1),'*r')
% % figure(4);
% for i=1:12
%     hold on
%  plot(t,rho(:,10,i))
% hold off
% end
% figure(5)
% plot(1:12, rho_unif(:,10),'*r')
% figure(6)
% plot(1:12,Z.*rho_unif(:,10)+rho_unif(:,10),'*r')
% for i=1:12
% rho_avg=rho(:,10,i).*model_end
% end

subplot(3,1,3)%figure(7)
plot(1:12, zeta_s, '*')