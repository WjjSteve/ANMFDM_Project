close all;
clear;

m = 201;
r_star = 0.1;
x1 = -1;
x2 = 1;
x0 = 0;
y1 = -1;
y2 = 1;
y0 = 0;
h = (x2-x1)/(m-1);
CFL =1.44;
x = x1:h:x2;
y = y1:h:y2;
T = 0.5;
k = CFL*h;
E_x = zeros(m);
for i=1:m
    for j=1:m
        H_z(i,j) = exp(-((x(i)-x0)/r_star)^2-((y(j)-y0)/r_star)^2);
    end
end
E_y = zeros(m);

figure(1)
surf(x,y,E_x);
hold on
surf(x,y,H_z);
hold on
surf(x,y,E_y);
legend("u^{(1)}","u^{(2)}","u^{(3)}");

u1 = [];
u2 = [];
u3 = [];
for i=1:m
    u1 = [u1,E_x(i,:)];
    u2 = [u2,H_z(i,:)];
    u3 = [u3,E_y(i,:)];
end
v = [u1';u2';u3'];

A = [0 0 0;0 0 -1;0 -1 0];
B = [0 1 0;1 0 0;0 0 0];
C = [1 0 0;0 1 0;0 0 1];
SBP4;

Dy = sparse(kron(eye(m),D1));
Dx = sparse(kron(D1,eye(m)));
HyI = sparse(kron(eye(m),HI));
HxI = sparse(kron(HI,eye(m)));
H_lI = sparse(HxI*HyI);
H_l = H_lI^(-1);
e_E = sparse(kron(e_m,eye(m)));
e_W = sparse(kron(e_1,eye(m)));
e_S = sparse(kron(eye(m),e_1));
e_N = sparse(kron(eye(m),e_m));
L_W = sparse(kron([0 1 0],e_W'));
L_E = sparse(kron([0 0 1],e_E'));
L_S = sparse(kron([1 0 0],e_S'));
L_N = sparse(kron([0 1 0],e_N'));
L_W = L_W(1:end-1,:);
L = sparse([L_E;L_W;L_N;L_S]);
C_l = sparse(kron(C,speye(m*m)));
H_hatI = sparse(kron(eye(3),H_lI));
KI = sparse(H_hatI*C_l);
P = sparse(speye(3*m*m)-KI*L'*(L*KI*L')^(-1)*L);
A_r = sparse(P*C_l^(-1)*(kron(A,Dx)+kron(B,Dy))*P);
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

U = v(1:m*m);
V = v(2*m*m+1:3*m*m);
THETA = Dx*U+Dy*V;
norm = sqrt(THETA'*H_l*THETA);
E_x = [];
H_z = [];
E_y = [];
u1 = v(1:m*m)';
u2 = v(m*m+1:2*m*m)';
u3 = v(2*m*m+1:3*m*m)';
for i=1:m
    E_x = [E_x;u1((i-1)*m+1:i*m)];
    H_z = [H_z;u2((i-1)*m+1:i*m)];
    E_y = [E_y;u3((i-1)*m+1:i*m)];
end
figure(2)
surf(x,y,E_x);
hold on
surf(x,y,H_z);
hold on
surf(x,y,E_y);
legend("u^{(1)}","u^{(2)}","u^{(3)}");


%5.0955e-06  5.0955e-06 2.6137e-04 0.0031     0.0056     0.0148
%3.7296e-10  4.1684e-08 4.5510e-06 5.8737e-04 0.0018     0.0064
%3.5921e-14  1.6302e-09 2.0032e-06 1.2300e-04 3.8956e-04 0.0014

