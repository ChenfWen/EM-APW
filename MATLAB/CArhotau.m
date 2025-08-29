function D = CArhotau(X)

global i phi gammas psi r alpha beta m u_CArhotau total_N avA avrho ame eps ksi num avtau

  N = X; 
  
  w(i)=beta*(avA*1.2*(N^eps))^(1/beta)*(alpha/r)^(alpha/beta)*((1-alpha-beta)/m)^((1-alpha-beta)/beta);
  %tau=RC(i)/w(i);
  w2(i)=((1-avtau)^(1+psi))*w(i);
  m1=272*(w2(i)^phi)+(w(i)^phi);
  c=( (1+psi)^(-(1+psi)) )*(psi^psi);
  D(i)= total_N*m1*avrho*((c*gammas*exp(ame(i)/(num(i,2)^ksi)*(N^ksi)-u_CArhotau))^phi)-N;

   
D = [D(i)];
