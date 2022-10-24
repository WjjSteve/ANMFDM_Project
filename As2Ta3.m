%close all;
%clear;

m = 401;
r_star = 0.1;
x_l = -1;
x_r = 1;
x_I = (x_r+x_l)/2;
l = (x_r-x_l)/2;
h = l/(m-1);
x = x_l:h:x_r;
xl = x_l:h:x_I;
xr = x_I:h:x_r;

T = 1.4;
CFL = 2.04;
k = CFL*h;
theta1 = exp(-(x+0.5).^2/(r_star)^2);
theta2 = -exp(-(x+0.5).^2/(r_star)^2);

u1 = (theta2-theta1)';
u2 = (theta2+theta1)';
figure(1)
plot(x,u1);
hold on
plot(x,u2);
legend("u^{(1)}","u^{(2)}");
%saveas(gcf,"a2t3ini.png");

v = [u1(1:m);u2(1:m);u1(m:2*m-1);u2(m:2*m-1)];
A = [0 1;1 0];
mul = 4;
mur = 1;
etal = sqrt(mul);
etar = sqrt(mur);
Cl = [mul 0;0 1];
Cr = [mur 0;0 1];
SBP4;

A_hat = [A,zeros(2,2);zeros(2,2),A];
C = kron([Cl,zeros(2,2);zeros(2,2),Cr],eye(m));
L = [kron([1 0],e_1'),zeros(1,2*m);zeros(1,2*m),kron([1 0],e_m');kron(eye(2),e_m'),-kron(eye(2),e_1')];
KI = C\kron(eye(4),HI); 
P = eye(4*m)-KI*L'*(L*KI*L')^(-1)*L;
A_r = P*C^(-1)*kron(A_hat,D1)*P;
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

u1 = [v(1:m);v(2*m+2:3*m)];
u41 = u1;
u2 = [v(m+1:2*m);v(3*m+2:4*m)];
u42 = u2;
figure(2)
plot(x,u1);
hold on
plot(x,u2);
legend("u^{(1)}","u^{(2)}");
%saveas(gcf,"a2t3t14.png");

a = [];
b = [];
for i = 2:400
    j = (i-1)*2;
    a = [a;u41(j+1)];
    b = [b;u42(j+1)];
end
a = [u41(1);a;u41(801)];
b = [u42(1);b;u42(801)];
%0.0206  0.0023 2.7066e-04
sqrt(0.02)*norm([a;b]-[u11;u12]);