function D = Crhoametau(X)

global i phi gammas psi r alpha beta m u_Crhoametau total_N A avrho avame eps ksi avtau num

  N = X; 
  
  w(i)=beta*(A(i)*(N^eps)/(num(i,2)^eps))^(1/beta)*(alpha/r)^(alpha/beta)*((1-alpha-beta)/m)^((1-alpha-beta)/beta);
  w2(i)=((1-avtau)^(1+psi))*w(i);
  m1=272*(w2(i)^phi)+(w(i)^phi);
  c=( (1+psi)^(-(1+psi)) )*(psi^psi);
  D(i)= total_N*m1*avrho*((c*gammas*exp(avame*(N^ksi)-u_Crhoametau))^phi)-N;

   
D = [D(i)];
