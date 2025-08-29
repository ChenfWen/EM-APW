function D = Ctau(X)

global i phi gammas psi r alpha beta m u_Ctau total_N A rho ame eps ksi num avtau

  N = X;  
  
  w(i)=beta(i)*(A(i)*(N^eps)/(num(i,2)^eps))^(1/beta(i))*(alpha(i)/r)^(alpha(i)/beta(i))*((1-alpha(i)-beta(i))/m(i))^((1-alpha(i)-beta(i))/beta(i));
  w2(i)=((1-avtau)^(1+psi))*w(i);
  m1=272*(w2(i)^phi)+w(i)^phi;
  c=( (1+psi)^(-(1+psi)) )*(psi^psi);
  D(i)= total_N*m1*rho(i)*((gammas*c*exp(ame(i)/(num(i,2)^ksi)*(N^ksi)-u_Ctau))^phi)-N;
  
D = [D(i)];

