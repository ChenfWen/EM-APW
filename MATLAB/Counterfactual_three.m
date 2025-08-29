clear,clc
clear all
global i phi gammas psi r alpha beta m u total_N A rho ame eps ksi num mrr N_huji tau
num=xlsread('data_ctrl.xls');
phi=1/0.3510;
gammas=gamma(1-1/phi);
psi=1.52; 
r=0.20;
alpha=0.486;
beta=0.504;
m=92476.196*(1-alpha-beta)/exp(6.174);
eps=0.00;
ksi=0.00;

N=num(:,2);
rho=num(:,3);
A=num(:,4);
N_huji=num(:,5);
total_N=sum(N);
mrr=N_huji/total_N;
u=10;


for i=1:length(num(:,1))
  w(i,1)=beta*A(i)^(1/beta)*(alpha/r)^(alpha/beta)*((1-alpha-beta)/m)^((1-alpha-beta)/beta);
  c1=exp(u)/gammas;
  c2=( (1+psi)^(-(1+psi)) )*(psi^psi);
  ame(i,1)=log(((mrr(i)/rho(i))^(1/phi))*c1/c2/w(i,1));
  w2(i,1)=((N(i)/total_N/rho(i)*((c1/c2/exp(ame(i,1)))^phi)-w(i,1)^phi)/272)^(1/phi);
  tau(i,1)=1-(w2(i)/w(i))^(1/(1+psi));
  RC(i,1)=tau(i)*w(i);
end

 

% Counterfactual A rho tau
count = 0; 
step = 1;
totalN_CArhotau = 0;
global u_CArhotau avA avrho avtau
u_CArhotau =u;
while abs(total_N-totalN_CArhotau) > 10000

    for i = 1:length(num(:,1))
       
         avA=sum(A.*N/total_N);
         avrho=sum(rho.*N/total_N);
         avtau=sum(tau.*N/total_N);
 
       if totalN_CArhotau == 0
            X0 = N(i);
        else
            X0 = max(N_CArhotau(i),1);
        end
        
        F=@(R)CArhotau(R);
        options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off');
        [X,Eval] = fsolve(F,X0,options);
    
        N_CArhotau(i) = max(0,real(X));
    end

 
    totalN_CArhotau = sum(N_CArhotau);
    count = count + 1;
    if count > 10
        step = step/10;
        count = 0;
    end
    
    if total_N-totalN_CArhotau > 10000
        u_CArhotau = u_CArhotau -step;
    else
        u_CArhotau = u_CArhotau +step;
    end
    disp(total_N-totalN_CArhotau)
end

totalN_CArhotau;
maxN_CArhotau=max(N_CArhotau);
u_CArhotau;
alloArhotau=sum(abs(N_CArhotau'-N))/2/total_N;







% Counterfactual rho ame tau
count = 0; 
step = 1;
totalN_Crhoametau = 0;
global u_Crhoametau avame
u_Crhoametau =u;
while abs(total_N-totalN_Crhoametau) > 10000

    for i = 1:length(num(:,1))
       
       avame=sum(ame.*N/total_N);
            
       if totalN_Crhoametau == 0
            X0 = N(i);
        else
            X0 = max(N_Crhoametau(i),1);
        end
        
        F=@(R)Crhoametau(R);
        options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off');
        [X,Eval] = fsolve(F,X0,options);
    
        N_Crhoametau(i) = max(0,real(X));
    end

 
    totalN_Crhoametau = sum(N_Crhoametau);
    count = count + 1;
    if count > 10
        step = step/10;
        count = 0;
    end
    
    if total_N-totalN_Crhoametau > 10000
        u_Crhoametau = u_Crhoametau -step;
    else
        u_Crhoametau = u_Crhoametau +step;
    end
    disp(total_N-totalN_Crhoametau)
end

totalN_Crhoametau;
maxN_Crhoametau=max(N_Crhoametau);
u_Crhoametau;
allorhoametau=sum(abs(N_Crhoametau'-N))/2/total_N;






% Counterfactual A ame tau
count = 0; 
step = 1;
totalN_CAametau = 0;
global u_CAametau
u_CAametau =u;
while abs(total_N-totalN_CAametau) > 10000

    for i = 1:length(num(:,1))
       
        
           
       if totalN_CAametau == 0
            X0 = N(i);
        else
            X0 = max(N_CAametau(i),1);
        end
        
        F=@(R)CAametau(R);
        options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off');
        [X,Eval] = fsolve(F,X0,options);
    
        N_CAametau(i) = max(0,real(X));
    end

 
    totalN_CAametau = sum(N_CAametau);
    count = count + 1;
    if count > 10
        step = step/10;
        count = 0;
    end
    
    if total_N-totalN_CAametau > 10000
        u_CAametau = u_CAametau -step;
    else
        u_CAametau = u_CAametau +step;
    end
    disp(total_N-totalN_CAametau)
end

totalN_CAametau;
maxN_CAametau=max(N_CAametau);
u_CAametau;
alloAametau=sum(abs(N_CAametau'-N))/2/total_N;






% Counterfactual A ame rho
count = 0; 
step = 1;
totalN_CAamerho = 0;
global u_CAamerho
u_CAamerho =u;
while abs(total_N-totalN_CAamerho) > 10000

    for i = 1:length(num(:,1))
       
        
           
       if totalN_CAamerho == 0
            X0 = N(i);
        else
            X0 = max(N_CAamerho(i),1);
        end
        
        F=@(R)CAamerho(R);
        options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off');
        [X,Eval] = fsolve(F,X0,options);
    
        N_CAamerho(i) = max(0,real(X));
    end

 
    totalN_CAamerho = sum(N_CAamerho);
    count = count + 1;
    if count > 10
        step = step/10;
        count = 0;
    end
    
    if total_N-totalN_CAamerho > 10000
        u_CAamerho = u_CAamerho -step;
    else
        u_CAamerho = u_CAamerho +step;
    end
    disp(total_N-totalN_CAamerho)
end

totalN_CAamerho;
maxN_CAamerho=max(N_CAamerho);
u_CAamerho;
alloAamerho=sum(abs(N_CAamerho'-N))/2/total_N;





for i = 1:length(N)
    lprob(i) = log((length(N) - i + 1)/length(N));
end
N1=sort(log(N));
N2=sort(log(N_Crhoametau));
N3=sort(log(N_CAametau));
N4=sort(log(N_CArhotau));
N5=sort(log(N_CAamerho));






a=[92,158,173]./255;
b=[210,204,161]./255;
subplot(1,1,1)
hold on
set(gca,'YGrid','on')
plot(N1,lprob,'Color',a,'LineWidth',2,'Linestyle','-');
plot(N2,lprob,'Color',b,'LineWidth',2,'Linestyle','--');
legend('Actual','Only efficiency','Location','SouthWest');
axis([11 19 -6 0]) 
xlabel('ln(population)')
ylabel('ln(prob > population)')
title({['Panel A.Utility=' num2str(u_Crhoametau)];['Reallocation=' num2str(allorhoametau)]},'Fontsize',10)


subplot(1,1,1)
hold on
set(gca,'YGrid','on')
plot(N1,lprob,'Color',a,'LineWidth',2,'Linestyle','-');
plot(N3,lprob,'Color',b,'LineWidth',2,'Linestyle','--');
legend('Actual','Only APW','Location','SouthWest');
axis([11 19 -6 0]) 
xlabel('ln(population)')
ylabel('ln(prob > population)')
title({['Panel B.Utility=' num2str(u_CAametau)];['Reallocation=' num2str(alloAametau)]},'Fontsize',10)

subplot(1,1,1)
hold on
set(gca,'YGrid','on')
plot(N1,lprob,'Color',a,'LineWidth',2,'Linestyle','-');
plot(N4,lprob,'Color',b,'LineWidth',2,'Linestyle','--');
legend('Actual','Only amenities','Location','SouthWest');
axis([11 19 -6 0]) 
xlabel('ln(population)')
ylabel('ln(prob > population)')
title({['Panel C.Utility=' num2str(u_CArhotau)];['Reallocation=' num2str(alloArhotau)]},'Fontsize',10)

subplot(1,1,1)
hold on
set(gca,'YGrid','on')
plot(N1,lprob,'Color',a,'LineWidth',2,'Linestyle','-');
plot(N5,lprob,'Color',b,'LineWidth',2,'Linestyle','--');
legend('Actual','Only migration friction','Location','SouthWest');
axis([11 19 -6 0]) 
xlabel('ln(population)')
ylabel('ln(prob > population)')
title({['Panel D.Utility=' num2str(u_CAamerho)];['Reallocation=' num2str(alloAamerho)]},'Fontsize',10)
