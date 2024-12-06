function [lmp, lambdapart]=lambda(q1,q2,rho,Risk)

global mu N0 epsilon nst rst NST
X=zeros(NST,rst);
 i=1:NST;
    for j = 1:rst
        isj=i+(j-1)*NST;
     X(i,j)=Risk(isj);
    end

%Definition of the Denominator for the risk groups:
%Low risk group rst=1 and high risk group rst=2
Dpart=zeros(NST,rst);
D=zeros(1,rst);
for j=1:rst
    for i=1:NST
        Dpart(i,j)= rho(i,j)*X(i,j); %each elemet in the Denominator
    end
D(j)=sum(Dpart(:,j));    %the whole denominator for lambda
end

%The Denominator for the G-matrix
        DTot=sum(D(:)); 

%The G-matrix
for i=1:rst
    for j=1:rst
        if i==j
            Delta(i,j)=1;   %Defining Delta Function
        else Delta(i,j)=0;
        end
         G(i,j)= epsilon*Delta(i,j)+(1-epsilon)*D(j)/DTot;
    end
end
lambdapart1=zeros(NST,rst,rst);lambdapart2=zeros(NST,rst,rst); lambdaT=zeros(rst,rst); lambdaTot=zeros(1,rst);

for i=1:rst
       for j=1:rst
            for beta=2:nst
           lambdapart1(beta,i,j)= rho(1,i)*q1(beta,j,i)*G(i,j)*rho(beta,j)*X(beta,j)/D(1,j);%one element of the lambda expression
            end
             for gamma=(nst+2):NST
           lambdapart2(gamma,i,j)= rho(1,i)*q2(gamma,j,i)*G(i,j)*rho(gamma,j)*X(gamma,j)/D(1,j);%one element of the lambda expression
             end
         lambdaTS(i,j)=sum(lambdapart1(:,i,j));
         lambdaTB(i,j)=sum(lambdapart2(:,i,j));
     lambdaT(i,j)=sum(lambdapart1(:,i,j))+sum(lambdapart2(:,i,j));
        end
     lambdaTot(i)= sum(lambdaT(i,:));
 end
lambdapart=lambdapart1+lambdapart2;
lmp=lambdaTot;
% 
% lambdapart=zeros(nst,rst,rst); lambdaT=zeros(rst,rst); lambdaTot=zeros(1,rst);
%  
% for i=1:rst
%      for j=1:rst
%          for alpha=2:nst
%            lambdapart(alpha,i,j)= rho(1,i)*q(alpha,j,i)*G(i,j)*rho(alpha,j)*X(alpha,j)/D(j);%one element of the lambda expression
%          end
%      lambdaT(i,j)=sum(lambdapart(:,i,j));
%      end
%      lambdaTot(i)= sum(lambdaT(i,:));
% end
% lmp=lambdaTot;
% end
% 
