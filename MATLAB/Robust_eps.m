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



%     for k=1:3
%     global u_m mp
%     mp=unifrnd(0.01,300*k,[273,1]).*m;
%     count = 0; 
%     step = 1;
%     totalN_m = 0;
%     u_m =u;
%     while abs(total_N-totalN_m) > 10000
% 
%     for i = 1:length(num(:,1))
%         
%        if totalN_m == 0
%             X0 = N(i);
%         else
%             X0 = max(N_m(i),1);
%         end
%         
%         F=@(R)dm(R);
%         options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off');
%         [X,Eval] = fsolve(F,X0,options);
%     
%         N_m(i) = max(0,real(X));
%     end
% 
%  
%     totalN_m = sum(N_m);
%     count = count + 1;
%     if count > 10
%         step = step/10;
%         count = 0;
%     end
%     
%     if total_N-totalN_m > 10000
%         u_m = u_m -step;
%     else
%         u_m = u_m +step;
%     end
%     disp(total_N-totalN_m)
% end
% totalN_m;
% maxN_m(k)=max(N_m);
% u_mr(k)=u_m;
% N_mr(:,k)=(N_m)';
% end
% 
% for i = 1:length(N)
%     lprob(i) = log((length(N) - i + 1)/length(N));
% end
% N=sort(log(N));
% N1=sort(log(N_mr(:,1)));
% N2=sort(log(N_mr(:,2)));
% N3=sort(log(N_mr(:,3)));
% 
% subplot(1,1,1)
% hold on
% set(gca,'YGrid','on')
% plot(N,lprob,'-r','LineWidth',1);
% plot(N1,lprob,'--b','LineWidth',1);
% plot(N2,lprob,'--g','LineWidth',1);
% plot(N3,lprob,'--k','LineWidth',1);
% axis([12 18 -6 0]) 
% legend('m at 10% level, utility = 9.9495 ','m at 50% level, utility = 9.9770 ','Actual: m, utility = 10','Location','SouthWest');
% xlabel('ln(population)')
% ylabel('ln(prob > population)')


phir=[2,2.5,phi,3,4];


for k = 1:5
    global u_phi phip gammasr
    phip=phir(k);
    gammasr=gamma(1-1/phip);
    count = 0; 
    step = 1;
    totalN_phi = 0;
    u_phi =u;
    while abs(total_N-totalN_phi) > 10000

    for i = 1:length(num(:,1))
        
       if totalN_phi == 0
            X0 = N(i);
        else
            X0 = max(N_phi(i),1);
        end
        
        F=@(R)dphi(R);
        options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off');
        [X,Eval] = fsolve(F,X0,options);
    
        N_phi(i) = max(0,real(X));
    end

 
    totalN_phi = sum(N_phi);
    count = count + 1;
    if count > 10
        step = step/10;
        count = 0;
    end
    
    if total_N-totalN_phi > 10000
        u_phi = u_phi -step;
    else
        u_phi = u_phi +step;
    end
    disp(total_N-totalN_phi)
end
totalN_phi;
maxN_phi(k)=max(N_phi);
u_phir(k)=u_phi;
N_phir(:,k)=(N_phi)';
end

for i = 1:length(N)
    lprob(i) = log((length(N) - i + 1)/length(N));
end
N1=sort(log(N_phir(:,1)));
N2=sort(log(N_phir(:,2)));
N3=sort(log(N_phir(:,3)));
N4=sort(log(N_phir(:,4)));
N5=sort(log(N_phir(:,5)));

a=[92,158,173]./255;
b=[210,204,161]./255;
c=[206,190,190]./255;
d=[237,177,131]./255;
e=[239,111,108]./255;
subplot(1,1,1)
hold on
set(gca,'YGrid','on')
plot(N1,lprob,'Color',c,'LineWidth',1.5,'Linestyle','--');
plot(N2,lprob,'Color',b,'LineWidth',1.5,'Linestyle','--');
plot(N3,lprob,'Color',a,'LineWidth',1.5,'Linestyle','-');
plot(N4,lprob,'Color',d,'LineWidth',1.5,'Linestyle','--');
plot(N5,lprob,'Color',e,'LineWidth',1.5,'Linestyle','--');
axis([10 20 -6 0]) 
legend(['\epsilon = 2, utility =',num2str(u_phir(1))],['\epsilon = 2.5, utility =',num2str(u_phir(2))],'Actual: \epsilon, utility = 10',['\epsilon = 3, utility =',num2str(u_phir(4))],['\epsilon = 4, utility =',num2str(u_phir(5))],'Location','SouthWest');
xlabel('ln(population)')
ylabel('ln(prob > population)')





%     global u_theta theta
%     for j=1:3
%     theta=normrnd(0.6,j/30,275,1);
%     count = 0; 
%     step = 1;
%     totalN_theta = 0;
%     u_theta =u;
%     while abs(total_N-totalN_theta) > 10000
% 
%     for i = 1:length(num(:,1))
%         
%        if totalN_theta == 0
%             X0 = N(i);
%         else
%             X0 = max(N_theta(i),1);
%         end
%         
%         F=@(R)dtheta(R);
%         options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off');
%         [X,Eval] = fsolve(F,X0,options);
%     
%         N_theta(i) = max(0,real(X));
%     end
% 
%  
%     totalN_theta = sum(N_theta);
%     count = count + 1;
%     if count > 10
%         step = step/10;
%         count = 0;
%     end
%     
%     if total_N-totalN_theta > 10000
%         u_theta = u_theta -step;
%     else
%         u_theta = u_theta +step;
%     end
%     disp(total_N-totalN_theta)
% end
% totalN_theta;
% maxN_theta=max(N_theta);
% u_thetar(:,j)=u_theta;
% N_thetar(:,j)=(N_theta)';
% end
% for i = 1:length(N)
%     lprob(i) = log((length(N) - i + 1)/length(N));
% end
% N1=sort(log(N));
% N2=sort(log(N_thetar(:,1)));
% N3=sort(log(N_thetar(:,2)));
% N4=sort(log(N_thetar(:,3)));
% 
% subplot(1,1,1)
% hold on
% set(gca,'YGrid','on')
% plot(N1,lprob,'-r','LineWidth',1);
% plot(N2,lprob,'--b','LineWidth',1);
% plot(N3,lprob,'--g','LineWidth',1);
% plot(N4,lprob,'--k','LineWidth',1);
% axis([10 20 -6 0]) 
% legend('Actual:f, utility = 10 ','f~N(0.6,0.033^2), utility = 10.0099 ','f~N(0.6,0.066^2), utility = 10.0246','f~N(0.6,0.1^2), utility = 10.0544 ','Location','SouthWest');
% 
% xlabel('ln(population)')
% ylabel('ln(prob > population)')
% 
% 
