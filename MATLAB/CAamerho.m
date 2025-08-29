function D = CAamerho(X)

global i phi gammas psi r alpha beta m u_CAamerho total_N avA avrho avame eps ksi tau

  N = X; 
  
  w(i)=beta*(avA*(N^eps))^(1/beta)*(alpha/r)^(alpha/beta)*((1-alpha-beta)/m)^((1-alpha-beta)/beta);
  w2(i)=((1-tau(i))^(1+psi))*w(i);
  m1=272*(w2(i)^phi)+(w(i)^phi);
  c=( (1+psi)^(-(1+psi)) )*(psi^psi);
  D(i)= total_N*m1*avrho*((c*gammas*exp(avame*(N^ksi)-u_CAamerho))^phi)-N;

   
D = [D(i)];
