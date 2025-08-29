clear,clc
clear all
global i phi gammas psi r alpha beta m u total_N A rho ame eps ksi num mrr N_huji tau
for k=1:5
num=xlsread('data_ctrl.xls');
pop=xlsread('pop_education.xls');
phi=1/0.3510;
gammas=gamma(1-1/phi);
psi=1.52; 
r=0.20;
alpha=0.486;
beta=0.504;
m=92476.196*(1-alpha-beta)/exp(6.174);
eps=0.04;
ksi=-0.04;
phirc=[1 1.5 2 phi 3];
phir=phirc(k);
N1=num(:,2);
N=pop(:,3);
rho=num(:,3).^phir;
A=num(:,4);
N_huji=num(:,5);
total_N1=sum(N1);
total_N=sum(N);
mrr=N_huji/total_N1;
u=10;


for i=1:length(num(:,1))
  w(i,1)=beta*A(i)^(1/beta)*(alpha/r)^(alpha/beta)*((1-alpha-beta)/m)^((1-alpha-beta)/beta);
  c1=exp(u)/gammas;
  c2=( (1+psi)^(-(1+psi)) )*(psi^psi);
  ame(i,1)=log(((mrr(i)/(rho(i)^(1/phir)))^(1/phi))*c1/c2/w(i,1));
  w2(i,1)=((N(i)/total_N/(rho(i))*((c1/c2/exp(ame(i,1)))^phi)-w(i,1)^phi)/272)^(1/phi);
  tau(i,1)=1-(w2(i)/w(i))^(1/(1+psi));
  RC(i,1)=tau(i)*w(i);
end

num(:,2)=N;

% Counterfactual A
count = 0; 
step = 1;
totalN_CA = 0;
global u_CA avA
u_CA =u;
while abs(total_N-totalN_CA) > 10000

    for i = 1:length(num(:,1))
       
        
            avA=sum(A./(N1.^eps).*(N1/total_N1));
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

N_CAr(:,k)=N_CA;
u_CAr(k)=u_CA;
alloAr(k)=sum(abs(N_CA'-N))/2/total_N;







% Counterfactual rho
count = 0; 
step = 1;
totalN_Crho = 0;
global u_Crho avrho
u_Crho =u;

if ksi==-0.04
    u_Crho =12;
end

while abs(total_N-totalN_Crho) > 10000

    for i = 1:length(num(:,1))
       
        
            avrho=(sum(rho.^(1/phir).*N1/total_N1))^phir;   
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

N_Crhor(:,k)=N_Crho;
u_Crhor(k)=u_Crho;
allorhor(k)=sum(abs(N_Crho'-N))/2/total_N;







% Counterfactual ame
count = 0; 
step = 1;
totalN_Came = 0;
global u_Came avame
u_Came =u;
while abs(total_N-totalN_Came) > 10000

    for i = 1:length(num(:,1))
       
        
            avame=sum(ame./(N1.^ksi).*(N1/total_N1));   
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

N_Camer(:,k)=N_Came;
u_Camer(k)=u_Came;
alloamer(k)=sum(abs(N_Came'-N))/2/total_N;



% Counterfactual tau
count = 0; 
step = 1;
totalN_Ctau = 0;
global u_Ctau avtau
u_Ctau =u;
while abs(total_N-totalN_Ctau) > 10000

    for i = 1:length(num(:,1))
       
        
            avtau=sum(tau.*N1/total_N1);   
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

N_Ctaur(:,k)=N_Ctau;
u_Ctaur(k)=u_Ctau;
allotaur(k)=sum(abs(N_Ctau'-N))/2/total_N;






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
end


for i = 1:length(N)
    lprob(i) = log((length(N) - i + 1)/length(N));
end
N0=sort(log(N));
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
plot(N0,lprob,'Color',a,'LineWidth',2,'Linestyle','-' );
plot(N21,lprob,'Color',b,'LineWidth',2,'Linestyle','--');
legend('Actual','Avg. efficiency','Location','SouthWest');
axis([8 17 -6 0]) 
xlabel('ln(population)')
ylabel('ln(prob > population)')
title({['Panel A.Utility=' num2str(u_CA)];['Reallocation=' num2str(alloAr(1))]},'Fontsize',10)


subplot(1,1,1)
hold on
set(gca,'YGrid','on')
plot(N0,lprob,'Color',a,'LineWidth',1.5,'Linestyle','-');
plot(N31,lprob,'Color',b,'LineWidth',1.5,'Linestyle','--');
plot(N32,lprob,'Color',c,'LineWidth',1.5,'Linestyle','--');
plot(N33,lprob,'Color',d,'LineWidth',1.5,'Linestyle','--');
plot(N34,lprob,'Color',e,'LineWidth',1.5,'Linestyle','--');
plot(N35,lprob,'Color',f,'LineWidth',1.5,'Linestyle','--');
legend('Actual,u=10',['\chi = 1, u =',num2str(u_Crhor(1))],['\chi = 1.5, u =',num2str(u_Crhor(2))],['\chi = 2, u =',num2str(u_Crhor(3))],['\chi = \epsilon, u =',num2str(u_Crhor(4))],['\chi = 3, u =',num2str(u_Crhor(5))],'Location','SouthWest');
axis([5 19 -6 0]) 
xlabel('ln(population)')
ylabel('ln(prob > population)')
title('Panel B.Avg.APW','Fontsize',10)


subplot(1,1,1)
hold on
set(gca,'YGrid','on')
plot(N0,lprob,'Color',a,'LineWidth',2,'Linestyle','-' );
plot(N41,lprob,'Color',b,'LineWidth',2,'Linestyle','--');
legend('Actual','Avg. amenities','Location','SouthWest');
axis([8 17 -6 0]) 
xlabel('ln(population)')
ylabel('ln(prob > population)')
title({['Panel C.Utility=' num2str(u_Came)];['Reallocation=' num2str(alloamer(1))]},'Fontsize',10)


subplot(1,1,1)
hold on
set(gca,'YGrid','on')
plot(N0,lprob,'Color',a,'LineWidth',1.5,'Linestyle','-');
plot(N51,lprob,'Color',b,'LineWidth',1.5,'Linestyle','--');
plot(N52,lprob,'Color',c,'LineWidth',1.5,'Linestyle','--');
plot(N53,lprob,'Color',d,'LineWidth',1.5,'Linestyle','--');
plot(N54,lprob,'Color',e,'LineWidth',1.5,'Linestyle','--');
plot(N55,lprob,'Color',f,'LineWidth',1.5,'Linestyle','--');
legend('Actual,u=10',['\chi = 1, u =',num2str(u_Ctaur(1))],['\chi = 1.5, u =',num2str(u_Ctaur(2))],['\chi = 2, u =',num2str(u_Ctaur(3))],['\chi = \epsilon, u =',num2str(u_Ctaur(4))],['\chi = 3, u =',num2str(u_Ctaur(5))],'Location','SouthWest');
axis([5 17 -6 0]) 
xlabel('ln(population)')
ylabel('ln(prob > population)')
title('Panel D.Avg.migration friction','Fontsize',10)


