close all;
clear;

m = 201;
r_star = 0.1;
x_l = -1;
x_r = 1;
l = x_r-x_l;
h = l/(m-1);
x = x_l:h:x_r;

T = 0;
CFL = 1.77;
k = CFL*h;
theta1 = exp(-x.^2/(r_star)^2);
theta2 = -exp(-x.^2/(r_star)^2);
v=[(theta2-theta1)';(theta2+theta1)'];

A = [0 1;1 0];
SBP6;
L = [kron([1 0],e_1');kron([1 0],e_m')];
HII = kron(eye(size(A,1)),HI); 
P = eye(2*m)-HII*L'*(L*HII*L')^(-1)*L;
A_r = P*kron(A,D1)*P;
t = 0;

while t<T
    if t+k>T
        k=T-t;
    end
    w1=A_r*v;
    w2=A_r*(v+k/2*w1);
    w3=A_r*(v+k/2*w2);
    w4=A_r*(v+k*w3);
    v=v+k/6*(w1+2*w2+2*w3+w4);
    t=t+k;
end

figure(1)
plot(x,v(1:m));
hold on
plot(x,v(m+1:end));
legend("u^{(1)}","u^{(2)}");
%saveas(gcf,"t11.png");

thetaf1 = exp(-((x-(l-T))/r_star).^2);
thetaf2 = -exp(-((x+(l-T))/r_star).^2);
u = [(thetaf1-thetaf2)';(thetaf1+thetaf2)'];
figure(2)
plot(x,u(1:m));
hold on
plot(x,u(m+1:end));
e = u-v;

NOE = sqrt(h)*norm(e)
% L2E2 = [0.4022 0.2467 0.0876 0.0215]
% L2E4 = [0.2226 0.0625 0.0067 4.7231e-04]
% L2E6 = [0.1816 0.0382 0.0032 2.0047e-04]