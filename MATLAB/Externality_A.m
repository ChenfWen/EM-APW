clear,clc
clear all
global i phi gammas psi r alpha beta m u total_N A rho ame eps ksi num mrr N_huji tau
num=xlsread('data_ctrl.xls');
phi=1/0.3510;
gammas=gamma(1-1/phi);
psi=1.52; 
r=0.20;
alpha=0.485;
beta=0.504;
m=115813.26*(1-alpha-beta)/exp(6.177);
epsp=[0.02 0.04 0.06 0.08 0.10];
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



 for j = 1:5
     eps=epsp(j);
% Counterfactual A
count = 0; 
step = 0.1;
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
N_CAr(:,j)=N_CA;
u_CAr(j)=u_CA;







% Counterfactual rho
count = 0; 
step = 0.1;
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

N_Crhor(:,j)=N_Crho;
u_Crhor(j)=u_Crho;





% Counterfactual ame
count = 0; 
step = 0.1;
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
N_Camer(:,j)=N_Came;
u_Camer(j)=u_Came;




count = 0; 
step = 0.1;
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
N_Ctaur(:,j)=N_Ctau;
u_Ctaur(j)=u_Ctau;


end

for i = 1:length(N)
    lprob(i) = log((length(N) - i + 1)/length(N));
end
N1=sort(log(N));
N21=sort(log(N_CAr(:,1)));
N22=sort(log(N_CAr(:,2)));
N23=sort(log(N_CAr(:,3)));
N24=sort(log(N_CAr(:,4)));
N25=sort(log(N_CAr(:,5)));
N31=sort(log(N_Crhor(:,1)));
N32=sort(log(N_Crhor(:,2)));
N33=sort(log(N_Crhor(:,3)));
N34=sort(log(N_Crhor(:,4)));
N35=sort(log(N_Crhor(:,5)));
N41=sort(log(N_Camer(:,1)));
N42=sort(log(N_Camer(:,2)));
N43=sort(log(N_Camer(:,3)));
N44=sort(log(N_Camer(:,4)));
N45=sort(log(N_Camer(:,5)));
N51=sort(log(N_Ctaur(:,1)));
N52=sort(log(N_Ctaur(:,2)));
N53=sort(log(N_Ctaur(:,3)));
N54=sort(log(N_Ctaur(:,4)));
N55=sort(log(N_Ctaur(:,5)));




a=[92,158,173]./255;
b=[210,204,161]./255;
c=[206,190,190]./255;
d=[237,177,131]./255;
e=[239,111,108]./255;
f=[48,151,164]./255;


subplot(1,1,1)
hold on
set(gca,'YGrid','on')
plot(N1,lprob,'Color',a,'LineWidth',1.5,'Linestyle','-');
plot(N21,lprob,'Color',b,'LineWidth',1.5,'Linestyle','--');
plot(N22,lprob,'Color',c,'LineWidth',1.5,'Linestyle','--');
plot(N23,lprob,'Color',d,'LineWidth',1.5,'Linestyle','--');
plot(N24,lprob,'Color',e,'LineWidth',1.5,'Linestyle','--');
plot(N25,lprob,'Color',f,'LineWidth',1.5,'Linestyle','--');
legend('Actual,u=10',['\delta = 0.02, u =',num2str(u_CAr(1))],['\delta = 0.04, u =',num2str(u_CAr(2))],['\delta = 0.06, u =',num2str(u_CAr(3))],['\delta = 0.08, u =',num2str(u_CAr(4))],['\delta = 0.10, u =',num2str(u_CAr(5))],'Location','SouthWest');
axis([11 19 -6 0]) 
xlabel('ln(population)')
ylabel('ln(prob > population)')
title(['Panel A.Avg.efficiency'],'Fontsize',10)


subplot(1,1,1)
hold on
set(gca,'YGrid','on')
plot(N1,lprob,'Color',a,'LineWidth',1.5,'Linestyle','-');
plot(N31,lprob,'Color',b,'LineWidth',1.5,'Linestyle','--');
plot(N32,lprob,'Color',c,'LineWidth',1.5,'Linestyle','--');
plot(N33,lprob,'Color',d,'LineWidth',1.5,'Linestyle','--');
plot(N34,lprob,'Color',e,'LineWidth',1.5,'Linestyle','--');
plot(N35,lprob,'Color',f,'LineWidth',1.5,'Linestyle','--');
legend('Actual,u=10',['\delta = 0.02, u =',num2str(u_Crhor(1))],['\delta = 0.04, u =',num2str(u_Crhor(2))],['\delta = 0.06, u =',num2str(u_Crhor(3))],['\delta = 0.08, u =',num2str(u_Crhor(4))],['\delta = 0.10, u =',num2str(u_Crhor(5))],'Location','SouthWest');
axis([11 19 -6 0]) 
xlabel('ln(population)')
ylabel('ln(prob > population)')
title(['Panel B.Avg.APW'],'Fontsize',10)


subplot(1,1,1)
hold on
set(gca,'YGrid','on')
plot(N1,lprob,'Color',a,'LineWidth',1.5,'Linestyle','-');
plot(N41,lprob,'Color',b,'LineWidth',1.5,'Linestyle','--');
plot(N42,lprob,'Color',c,'LineWidth',1.5,'Linestyle','--');
plot(N43,lprob,'Color',d,'LineWidth',1.5,'Linestyle','--');
plot(N44,lprob,'Color',e,'LineWidth',1.5,'Linestyle','--');
plot(N45,lprob,'Color',f,'LineWidth',1.5,'Linestyle','--');
legend('Actual,u=10',['\delta = 0.02, u =',num2str(u_Camer(1))],['\delta = 0.04, u =',num2str(u_Camer(2))],['\delta = 0.06, u =',num2str(u_Camer(3))],['\delta = 0.08, u =',num2str(u_Camer(4))],['\delta = 0.10, u =',num2str(u_Camer(5))],'Location','SouthWest');
axis([11 19 -6 0]) 
xlabel('ln(population)')
ylabel('ln(prob > population)')
title(['Panel C.Avg.amenities'],'Fontsize',10)

subplot(1,1,1)
hold on
set(gca,'YGrid','on')
plot(N1,lprob,'Color',a,'LineWidth',1.5,'Linestyle','-');
plot(N51,lprob,'Color',b,'LineWidth',1.5,'Linestyle','--');
plot(N52,lprob,'Color',c,'LineWidth',1.5,'Linestyle','--');
plot(N53,lprob,'Color',d,'LineWidth',1.5,'Linestyle','--');
plot(N54,lprob,'Color',e,'LineWidth',1.5,'Linestyle','--');
plot(N55,lprob,'Color',f,'LineWidth',1.5,'Linestyle','--');
legend('Actual,u=10',['\delta = 0.02, u =',num2str(u_Ctaur(1))],['\delta = 0.04, u =',num2str(u_Ctaur(2))],['\delta = 0.06, u =',num2str(u_Ctaur(3))],['\delta = 0.08, u =',num2str(u_Ctaur(4))],['\delta = 0.10, u =',num2str(u_Ctaur(5))],'Location','SouthWest');
axis([11 19 -6 0]) 
xlabel('ln(population)')
ylabel('ln(prob > population)')
title(['Panel D.Avg.migration friction'],'Fontsize',10)


