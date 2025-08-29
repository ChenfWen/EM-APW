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

ekrs=[0.06 0.02 0.04 0.06];
ekri=[-0.06 0.02 0.04 0.06];

%%
    global u_ek epsp ksip avA avrho avame avtau
    epsp=ekrs(1);
    ksip=ekri(1);
    count = 0; 
    step = 0.10;
    totalN_ek = 0;
    u_ek =u;
    while abs(total_N-totalN_ek) > 10000
        avA=sum(A./(N.^epsp).*(N/total_N));
        avrho=sum(rho.*N/total_N);
        avame=sum(ame./(N.^ksip).*(N/total_N));
        avtau=sum(tau.*N/total_N); 

    for i = 1:length(num(:,1))
        
       if totalN_ek == 0
            X0 = N(i);
        else
            X0 = max(N_ek(i),1);
        end
        
        F=@(R)dek_no(R);
        options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off');
        [X,Eval] = fsolve(F,X0,options);
    
        N_ek(i) = max(0,real(X));
    end

 
    totalN_ek = sum(N_ek);
    count = count + 1;
    if count > 10
        step = step/10;
        count = 0;
    end
    
    if total_N-totalN_ek > 10000
        u_ek = u_ek -step;
    else
        u_ek = u_ek +step;
    end
    disp(total_N-totalN_ek)
    end
u_ekr(1)=u_ek;
N_ekr(:,1)=(N_ek)';

%%
   global u_ek epsp ksip avA avrho avame avtau
    epsp=ekrs(2);
    ksip=ekri(2);
    count = 0; 
    step = 0.1;
    totalN_ek = 0;
    u_ek =10.1;
    while abs(total_N-totalN_ek) > 10000
        avA=sum(A./(N.^epsp).*(N/total_N));
        avrho=sum(rho.*N/total_N);
        avame=sum(ame./(N.^ksip).*(N/total_N));
        avtau=sum(tau.*N/total_N); 


    for i = 1:length(num(:,1))
        
       if totalN_ek == 0
            X0 = N(i);
        else
            X0 = max(N_ek(i),1);
        end
        
        F=@(R)dek_no(R);
        options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off');
        [X,Eval] = fsolve(F,X0,options);
    
        N_ek(i) = max(0,real(X));
    end

 
    totalN_ek = sum(N_ek);
    count = count + 1;
    if count > 10
        step = step/10;
        count = 0;
    end
    
    if total_N-totalN_ek > 10000
        u_ek = u_ek -step;
    else
        u_ek = u_ek +step;
    end
    disp(total_N-totalN_ek)
    end
u_ekr(2)=u_ek;
N_ekr(:,2)=(N_ek)';


%%   
  global u_ek epsp ksip avA avrho avame avtau
    epsp=ekrs(3);
    ksip=ekri(3);
    count = 0; 
    step = 0.02;
    totalN_ek = 0;
    u_ek =10.33;
    while abs(total_N-totalN_ek) > 10000
        avA=sum(A./(N.^epsp).*(N/total_N));
        avrho=sum(rho.*N/total_N);
        avame=sum(ame./(N.^ksip).*(N/total_N));
        avtau=sum(tau.*N/total_N); 


    for i = 1:length(num(:,1))
        
       if totalN_ek == 0
            X0 = N(i);
        else
            X0 = max(N_ek(i),1);
        end
        
        F=@(R)dek_no(R);
        options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off');
        [X,Eval] = fsolve(F,X0,options);
    
        N_ek(i) = max(0,real(X));
    end

 
    totalN_ek = sum(N_ek);
    count = count + 1;
    if count > 10
        step = step/10;
        count = 0;
    end
    
    if total_N-totalN_ek > 10000
        u_ek = u_ek +step;
    else
        u_ek = u_ek -step;
    end
    disp(total_N-totalN_ek)
    end
u_ekr(3)=u_ek;
N_ekr(:,3)=(N_ek)';


%%   
    global u_ek epsp ksip avA avrho avame avtau
    epsp=ekrs(4);
    ksip=ekri(4);
    count = 0; 
    step = 0.1;
    totalN_ek = 0;
    u_ek =10.3;
    while abs(total_N-totalN_ek) > 10000
        avA=sum(A./(N.^epsp).*(N/total_N));
        avrho=sum(rho.*N/total_N);
        avame=sum(ame./(N.^ksip).*(N/total_N));
        avtau=sum(tau.*N/total_N); 


    for i = 1:length(num(:,1))
        
       if totalN_ek == 0
            X0 = N(i);
        else
            X0 = max(N_ek(i),1);
        end
        
        F=@(R)dek_no(R);
        options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off');
        [X,Eval] = fsolve(F,X0,options);
    
        N_ek(i) = max(0,real(X));
    end

 
    totalN_ek = sum(N_ek);
    count = count + 1;
    if count > 10
        step = step/10;
        count = 0;
    end
    
    if total_N-totalN_ek > 10000
        u_ek = u_ek +step;
    else
        u_ek = u_ek -step;
    end
    disp(total_N-totalN_ek)
    end
u_ekr(4)=u_ek;
N_ekr(:,4)=(N_ek)';


%%
for i = 1:length(N)
    lprob(i) = log((length(N) - i + 1)/length(N));
end
N1=sort(log(N_ekr(:,1)));
N2=sort(log(N_ekr(:,2)+1));
N3=sort(log(N_ekr(:,3)+1));
N4=sort(log(N_ekr(:,4)+1));
N5=sort(log(N));

a=[92,158,173]./255;
b=[210,204,161]./255;
c=[206,190,190]./255;
d=[237,177,131]./255;
e=[239,111,108]./255;
subplot(1,1,1)
hold on
set(gca,'YGrid','on')
plot(N5,lprob,'Color',a,'LineWidth',1.,'Linestyle','-');
plot(N1,lprob,'Color',b,'LineWidth',1.,'Linestyle','--','Marker','.','MarkerSize',8);
plot(N2,lprob,'Color',c,'LineWidth',1.,'Linestyle','--','Marker','.','MarkerSize',10);
plot(N3,lprob,'Color',d,'LineWidth',1.,'Linestyle','--','Marker','.','MarkerSize',12);
plot(N4,lprob,'Color',e,'LineWidth',1.,'Linestyle','--','Marker','.','MarkerSize',15);
axis([10 20 -6 0]) 
legend('Actual , utility = 10',['\delta = -\omega = 0.06, utility =',num2str(u_ekr(1))],['\delta = \omega = 0.02, utility =',num2str(u_ekr(2))],['\delta = \omega = 0.04, utility =',num2str(u_ekr(3))],['\delta = \omega = 0.06, utility =',num2str(u_ekr(4))], 'Location' ,'SouthWest');
xlabel('ln(population)')
ylabel('ln(prob > population)')
