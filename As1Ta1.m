m =201;
h = 1/(m-1);

A = [0 1;1 0];
SBP6;

%L = [kron([1 0],e_1');kron([0 1],e_m')];
L = [kron([1 -1],e_1');kron([1 1],e_m')];
HII = kron(eye(size(A,1)),HI); 
P = eye(2*m)-HII*L'*(L*HII*L')^(-1)*L;
A_r = P*kron(A,D1)*P;
e = eig(h*A_r);

figure(1)
plot(real(e),imag(e),'r*') %   Plot real and imaginary parts
xlabel('Real')
ylabel('Imaginary')
r = max(abs(e))
%saveas(gcf,"6SBPm51Cha.png");
