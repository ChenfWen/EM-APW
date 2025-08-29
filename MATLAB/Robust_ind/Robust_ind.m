clear,clc
clear all
global i phi gammas psi r alpha beta m u total_N A rho ame eps ksi num mrr N_huji tau
num=xlsread('data_ctrl.xls');
phi=1/0.3510;
gammas=gamma(1-1/phi);
psi=1.52; 
r=0.20;
eps=0.00;
ksi=0.00;

N=num(:,2);
rho=num(:,3);
A=num(:,4);
N_huji=num(:,5);
alpha=num(:,6);
beta=num(:,7);
rho=rho.*(1-alpha-beta)/(0.01);
m=92476.196*(1-alpha-beta)/exp(6.174);
total_N=sum(N);
mrr=N_huji/total_N;
u=10;


for i=1:length(num(:,1))
  w(i,1)=beta(i)*A(i)^(1/beta(i))*(alpha(i)/r)^(alpha(i)/beta(i))*((1-alpha(i)-beta(i))/m(i))^((1-alpha(i)-beta(i))/beta(i));
  c1=exp(u)/gammas;
  c2=( (1+psi)^(-(1+psi)) )*(psi^psi);
  ame(i,1)=log(((mrr(i)/rho(i))^(1/phi))*c1/c2/w(i,1));
  w2(i,1)=((N(i)/total_N/rho(i)*((c1/c2/exp(ame(i,1)))^phi)-w(i,1)^phi)/272)^(1/phi);
  tau(i,1)=1-(w2(i)/w(i))^(1/(1+psi));
  RC(i,1)=tau(i)*w(i);
end


% Counterfactual A
count = 0; 
step = 1;
totalN_CA = 0;
global u_CA avA
u_CA =u;
while abs(total_N-totalN_CA) > 10000

    for i = 1:length(num(:,1))
       
        
            avA=sum(A./(N.^eps).*(N/total_N));
       if totalN_CA == 0
            X0 = N(i);
        else
            X0 = max(N_CA(i),1);
        end
        
        F=@(R)CA(R);
        options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off');
        [X,Eval] = fsolve(F,X0,options);
    
        N_CA(i) = max(0,real(X));
    end

 
    totalN_CA = sum(N_CA);
    count = count + 1;
    if count > 10
        step = step/10;
        count = 0;
    end
    
    if total_N-totalN_CA > 10000
        u_CA = u_CA -step;
    else
        u_CA = u_CA +step;
    end
    disp(total_N-totalN_CA)
end

totalN_CA;
maxN_CA=max(N_CA);
u_CA;
alloA=sum(abs(N_CA'-N))/2/total_N;







% Counterfactual rho
count = 0; 
step = 1;
totalN_Crho = 0;
global u_Crho avrho
u_Crho =u;
while abs(total_N-totalN_Crho) > 10000

    for i = 1:length(num(:,1))
       
        
            avrho=sum(rho.*N/total_N);   
       if totalN_Crho == 0
            X0 = N(i);
        else
            X0 = max(N_Crho(i),1);
        end
        
        F=@(R)Crho(R);
        options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off');
        [X,Eval] = fsolve(F,X0,options);
        N_Crho(i) = max(0,real(X));
    end

 
    totalN_Crho = sum(N_Crho);
    count = count + 1;
    if count > 10
        step = step/10;
        count = 0;
    end
    
    if total_N-totalN_Crho > 10000
        u_Crho = u_Crho -step;
    else
        u_Crho = u_Crho +step;
    end
    disp(total_N-totalN_Crho)
end

totalN_Crho;
maxN_Crho=max(N_Crho);
u_Crho;
allorho=sum(abs(N_Crho'-N))/2/total_N;






% Counterfactual ame
count = 0; 
step = 1;
totalN_Came = 0;
global u_Came avame
u_Came =u;
while abs(total_N-totalN_Came) > 10000

    for i = 1:length(num(:,1))
       
        
            avame=sum(ame./(N.^ksi).*(N/total_N));   
       if totalN_Came == 0
            X0 = N(i);
        else
            X0 = max(N_Came(i),1);
        end
        
        F=@(R)Came(R);
        options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off');
        [X,Eval] = fsolve(F,X0,options);
    
        N_Came(i) = max(0,real(X));
    end

 
    totalN_Came = sum(N_Came);
    count = count + 1;
    if count > 10
        step = step/10;
        count = 0;
    end
    
    if total_N-totalN_Came > 10000
        u_Came = u_Came -step;
    else
        u_Came = u_Came +step;
    end
    disp(total_N-totalN_Came)
end

totalN_Came;
maxN_Came=max(N_Came);
u_Came;
alloame=sum(abs(N_Came'-N))/2/total_N;


% Counterfactual tau
count = 0; 
step = 1;
totalN_Ctau = 0;
global u_Ctau avtau
u_Ctau =u;
while abs(total_N-totalN_Ctau) > 10000

    for i = 1:length(num(:,1))
       
        
            avtau=sum(tau.*N/total_N);   
       if totalN_Ctau == 0
            X0 = N(i);
        else
            X0 = max(N_Ctau(i),1);
        end
        
        F=@(R)Ctau(R);
        options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off');
        [X,Eval] = fsolve(F,X0,options);
        N_Ctau(i) = max(0,real(X));
    end

 
    totalN_Ctau = sum(N_Ctau);
    count = count + 1;
    if count > 10
        step = step/10;
        count = 0;
    end
    
    if total_N-totalN_Ctau > 10000
        u_Ctau = u_Ctau -step;
    else
        u_Ctau = u_Ctau +step;
    end
    disp(total_N-totalN_Ctau)
end

totalN_Ctau;
maxN_Ctau=max(N_Ctau);
u_Ctau;
allotau=sum(abs(N_Ctau'-N))/2/total_N;





% % Counterfactual RC
% count = 0; 
% step = 1;
% totalN_CRC = 0;
% global u_CRC avRC
% u_CRC =u;
% while abs(total_N-totalN_CRC) > 10000
% 
%     for i = 1:length(num(:,1))
%        
%         
%             avRC=sum(RC.*N/total_N);   
%        if totalN_CRC == 0
%             X0 = N(i);
%         else
%             X0 = max(N_CRC(i),1);
%         end
%         
%         F=@(R)CRC(R);
%         options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off');
%         [X,Eval] = fsolve(F,X0,options);
%         N_CRC(i) = max(0,real(X));
%     end
% 
%  
%     totalN_CRC = sum(N_CRC);
%     count = count + 1;
%     if count > 10
%         step = step/10;
%         count = 0;
%     end
%     
%     if total_N-totalN_CRC > 10000
%         u_CRC = u_CRC -step;
%     else
%         u_CRC = u_CRC +step;
%     end
%     disp(total_N-totalN_CRC)
% end
% 
% totalN_CRC;
% maxN_CRC=max(N_CRC);
% u_CRC;
% 
% 





for i = 1:length(N)
    lprob(i) = log((length(N) - i + 1)/length(N));
end
N1=sort(log(N));
N2=sort(log(N_CA));
N3=sort(log(N_Crho));
N4=sort(log(N_Came));
N5=sort(log(N_Ctau));






a=[92,158,173]./255;
b=[210,204,161]./255;
subplot(1,1,1)
hold on
set(gca,'YGrid','on')
plot(N1,lprob,'Color',a,'LineWidth',2,'Linestyle','-' );
plot(N2,lprob,'Color',b,'LineWidth',2,'Linestyle','--');
legend('Actual','Avg. efficiency','Location','SouthWest');
axis([11 19 -6 0]) 
xlabel('ln(population)')
ylabel('ln(prob > population)')
title({['Panel A.Utility=' num2str(u_CA)];['Reallocation=' num2str(alloA)]},'Fontsize',10)


subplot(1,1,1)
hold on
set(gca,'YGrid','on')
plot(N1,lprob,'Color',a,'LineWidth',2,'Linestyle','-');
plot(N3,lprob,'Color',b,'LineWidth',2,'Linestyle','--');
legend('Actual','Avg. APW','Location','SouthWest');
axis([11 19 -6 0]) 
xlabel('ln(population)')
ylabel('ln(prob > population)')
title({['Panel B.Utility=' num2str(u_Crho)];['Reallocation=' num2str(allorho)]},'Fontsize',10)

subplot(1,1,1)
hold on
set(gca,'YGrid','on')
plot(N1,lprob,'Color',a,'LineWidth',2,'Linestyle','-');
plot(N4,lprob,'Color',b,'LineWidth',2,'Linestyle','--');
legend('Actual','Avg. amenities','Location','SouthWest');
axis([11 19 -6 0]) 
xlabel('ln(population)')
ylabel('ln(prob > population)')
title({['Panel C.Utility=' num2str(u_Came)];['Reallocation=' num2str(alloame)]},'Fontsize',10)

subplot(1,1,1)
hold on
set(gca,'YGrid','on')
plot(N1,lprob,'Color',a,'LineWidth',2,'Linestyle','-');
plot(N5,lprob,'Color',b,'LineWidth',2,'Linestyle','--');
legend('Actual','Avg. migration friction','Location','SouthWest');
axis([11 19 -6 0]) 
xlabel('ln(population)')
ylabel('ln(prob > population)')
title({['Panel D.Utility=' num2str(u_Ctau)];['Reallocation=' num2str(allotau)]},'Fontsize',10)


% N_CA=(N_CA)';
% N_Crho=(N_Crho)';
% N_Came=(N_Came)';
% N_Ctau=(N_Ctau)';
% 
% 
% 
% 
% citygroup = importdata('citygroup.txt');
%  group1(:,1)=citygroup(1:29,1);
%  group2(:,1)=citygroup(1:21,2);
%  group3(:,1)=citygroup(:,3);
%  
% 
%  for i = 1:length(group1)
%    for j = 1:length(N_Crho)
%        if num(j,1)==group1(i,1)
%            group1(i,2)=N_Crho(j);
%            group1(i,3)=N(j);
%        end
%    end
%  end
%  
%  for i = 1:length(group2)
%    for j = 1:length(N_Crho)
%        if num(j,1)==group2(i,1)
%            group2(i,2)=N_Crho(j);
%            group2(i,3)=N(j);
%        end
%    end
%  end
%  
%   for i = 1:length(group3)
%    for j = 1:length(N_Crho)
%        if num(j,1)==group3(i,1)
%            group3(i,2)=N_Crho(j);
%            group3(i,3)=N(j);
%        end
%    end
%   end
%   
% 
%   for i = 1:length(group3)
%       if group3(i,2)==0
%           group3(i,:)=[];
%       end
%   end
%   
% 
%   for i = 1:length(group1)
%     lprob1(i) = log((length(group1) - i + 1)/length(group1));
% end
%   N5=sort(log(group1(:,3)));
%   N6=sort(log(group1(:,2)));
%   subplot(1,1,1)
% hold on
% set(gca,'YGrid','on')
% plot(N5,lprob1,'Color',a,'LineWidth',2,'Linestyle','-' );
% plot(N6,lprob1,'Color',b,'LineWidth',2,'Linestyle','--');
% legend('Actual','Avg. APW','Location','SouthWest');
% axis([12 19 -3.5 0]) 
% xlabel('ln(population)')
% ylabel('ln(prob > population)')
% title('Group 1')
% 
% for i = 1:length(group2)
%     lprob2(i) = log((length(group2) - i + 1)/length(group2));
% end
%   N7=sort(log(group2(:,3)));
%   N8=sort(log(group2(:,2)));
%   subplot(1,1,1)
% hold on
% set(gca,'YGrid','on')
% plot(N7,lprob2,'Color',a,'LineWidth',2,'Linestyle','-' );
% plot(N8,lprob2,'Color',b,'LineWidth',2,'Linestyle','--');
% legend('Actual',['Avg. APW'],'Location','SouthWest');
% axis([12 19 -3.5 0]) 
% xlabel('ln(population)')
% ylabel('ln(prob > population)')
% title('Group 2')
% 
% 
% for i = 1:length(group3)
%     lprob3(i) = log((length(group3) - i + 1)/length(group3));
% end
%   N9=sort(log(group3(:,3)));
%   N10=sort(log(group3(:,2)));
%   subplot(1,1,1)
% hold on
% set(gca,'YGrid','on')
% plot(N9,lprob3,'Color',a,'LineWidth',2,'Linestyle','-' );
% plot(N10,lprob3,'Color',b,'LineWidth',2,'Linestyle','--');
% legend('Actual',['Avg. APW'],'Location','SouthWest');
% axis([12 19 -4.5 0]) 
% xlabel('ln(population)')
% ylabel('ln(prob > population)')
% title('Group 3')
