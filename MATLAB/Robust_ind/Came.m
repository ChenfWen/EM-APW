function D = Came(X)

global i phi gammas psi r alpha beta m u_Came total_N A rho avame eps ksi num tau

  N = X; 
  
  w(i)=beta(i)*(A(i)*(N^eps)/(num(i,2)^eps))^(1/beta(i))*(alpha(i)/r)^(alpha(i)/beta(i))*((1-alpha(i)-beta(i))/m(i))^((1-alpha(i)-beta(i))/beta(i));
  w2(i)=((1-tau(i))^(1+psi))*w(i);
  m1=272*(w2(i)^phi)+w(i)^phi;
  c=( (1+psi)^(-(1+psi)) )*(psi^psi);
  D(i)= total_N*m1*rho(i)*((gammas*c*exp(avame*(N^ksi)-u_Came))^phi)-N;
  
D = [D(i)];

