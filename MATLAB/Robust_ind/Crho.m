function D = Crho(X)

global i phi gammas psi r alpha beta m u_Crho total_N A avrho ame eps ksi num tau

  N = X;  
  
  w(i)=beta(i)*(A(i)*(N^eps)/(num(i,2)^eps))^(1/beta(i))*(alpha(i)/r)^(alpha(i)/beta(i))*((1-alpha(i)-beta(i))/m(i))^((1-alpha(i)-beta(i))/beta(i));
  w2(i)=((1-tau(i))^(1+psi))*w(i);
  m1=272*(w2(i)^phi)+w(i)^phi;
  c=( (1+psi)^(-(1+psi)) )*(psi^psi);
  D(i)= total_N*m1*avrho*((gammas*c*exp(ame(i)/(num(i,2)^ksi)*(N^ksi)-u_Crho))^phi)-N;
  
D = [D(i)];
