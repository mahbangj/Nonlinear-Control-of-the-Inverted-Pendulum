clc;
clear all;
%% 1
syms m1 m2 m3 g r1 r2
syms x1 x2 x3 x4 x5 x6

m1 = 2;
m2 = 0.1;
m3 = 0.1;
r1 = 0.1;
r2 = 0.1;
g = 9.8;


%A_1 = [0 1 0 0;0 0 (-m2*g*sin(x3)*cos(x3))/(x3(m1+m2*sin(x3)^2)) (m2*r*x4*sin(x3))/(m1+m2*sin(x3)^2);0 0 0 1;0 0 ((m1+m2)*g*sin(x3))/(r*x3*(m1+m2*sin(x3)^2)) (-m2*r*x4*sin(x3)*cos(x3))/(r*(m1+m2*sin(x3)^2))]
%B_1 = [0 1/(m1+m2*sin(x3)^2) 0 -cos(x3)/(r*(m1+m2*sin(x3)^2))]'
a = 2+0.1;
b = cos(x3)*cos(x5);
c = sin(x3-x5);
d = 0.1*a*sin(x3)^2 + 2*0.1*c^2 + 2*0.1;
e = 2*0.1^2 + 2*0.1 + 2*0.1*0.1 + 2*2*0.1 + 2*0.1*cos(2*x5);
f = 0.1*a*sin(2*x3) - 2*0.1*sin(2*(x5-x3));
h = 0.1*0.1*sin(x3)*cos(x5) + 2*0.1*c;

Q=(0.1*9.8*a*sin(x3)*cos(x3))/(-d*x3);
W=0.1*0.1*a*x4*sin(x3)/d;
R=0.1*0.1*0.1*x6*(sin(x5-2*x3)-sin(x5))/(-2*d);
T=9.8*e*sin(x3)/(2*0.1*d*x3);
Y=f*x4/(-2*d);
U=2*0.1*9.8*b*sin(x5)/(-0.1*d*x5);
I=0.1*h*x6/(-0.1*d);
O=2*9.8*a*b*sin(x3)/(-0.1*d*x3);
P=2*0.1*a*c*x4/(0.1*d);
S=2*9.8*a*((cos(x3))^2)*sin(x5)/(0.1*d*x5);
D=2*0.1*x6*sin(2*(x5-x3))/(-2*d);
A_n = [0 1 0 0 0 0;0 0 Q W 0 R;0 0 0 1 0 0;0 0 T Y U I;0 0 0 0 0 1;0 0 O P S D]; 
%A_n = [0 1 0 0 0 0;0 0 (m2*g*a*sin(x3)*cos(x3))/(-d*x3) m2*r1*a*x4*sin(x3)/d 0 m2*m3*r2*x6*(sin(x5-2*x3)-sin(x5))/(-2*d);0 0 0 1 0 0;0 0 g*e*sin(x3)/(2*r1*d*x3) f*x4/(-2*d) m1*m3*g*b*sin(x5)/(-r1*d*x5) r2*h*x6/(-r1*d);0 0 0 0 0 1;0 0 m1*g*a*b*sin(x3)/(-r2*d*x3) m1*r1*a*c*x4/(r2*d) m1*g*a*((cos(x3))^2)*sin(x5)/(r2*d*x5) m1*m3*x6*sin(2*(x5-x3))/(-2*d)]  

B_n = [0 (m2+m3*c^2)/d 0 (m2*cos(x3)-m3*c*sin(x5))/(-r1*d) 0 (m2*cos(x3)-m3*c*sin(x5))/(-r1*d)]'
A_nn=[0 1 0 0 0 0 1 0 0 0 0 0;0 0 Q W 0 R 0 1 0 0 0 0;0 0 0 1 0 0 0 0 1 0 0 0;0 0 T Y U I 0 0 0 1 0 0;0 0 0 0 0 1 0 0 0 0 1 0;0 0 O P S D 0 0 0 0 0 1]
%% 



%syms t1 t2 t3 t4 real
%syms td1 td2 td3 td4 real
syms xd1 xd2 xd3 xd4 xd5 xd6 real
syms x1 x2 x3 x4 x5 x6 real
syms u real

%X = vpa(simplifyFraction(pinv(A_n)*B_n),3);
%X = pinv(A_n)*B_n; % for state representation
%  
%  % x_dot = F(x,u)
%  F = [xd1 xd2 xd3 xd4 xd5 xd6 X'];
%  st = [x1 x2 x3 x4 x5 x6 xd1 xd2 xd3 xd4 xd5 xd6];
%  As = jacobian(F,st);
%  A_l = subs(As,[st, u],zeros(1,13));
%  A_l = vpa(simplify(A_l),3);
%  
%  Bs = jacobian(F,u);
%  B_l = subs(Bs,[st, u],zeros(1,13));
%  B_l = vpa(simplify(B_l),3);
%  
%  C_l = [0 0 1 0 0 0;0 0 0 0 1 0];
%  D_l = [0;0];
%  
%  A_l = double(A_l);
%  B_l = double(B_l);
%  
%  sys = ss(A_l,B_l,C_l,D_l);
%  
% 
%  x0 = [0 0 pi/6 0 0 0]';
%  initial(sys,x0,10)

%% 2
syms land
syms v1 v2 v3 v4 v5 v6
V =[v1;v2;v3;v4;v5;v6];
m1 = 2;
m2 = 0.1;
m3 = 0.1;
r1 = 0.1;
r2 = 0.1;
g = 9.8;
M = 2*m2^2 + m1*m3 + 2*m2*m3 + 2*m1*m2 + m1*m3;
A = [0 1 0 0 0 0; 0 0 (-((m2+m3)*g)/m1) 0 0 0; 0 0 0 1 0 0; 0 0 ((M*g)/(2*r1*m1*m2)) 0 (-(m3*g)/(r1*m2)) 0;
     0 0 0 0 0 1; 0 0 (-((m2+m3)*g)/(r2*m2)) 0 (((m2+m3)*g)/(r2*m2)) 0];
B = [0 ;(1/m1);0;(-1/(r1*m1));0;0];
%eig = eig(A)
%[V,D] = eig(A)
% det = det(land * eye(6) - A)
% solve(det,land)
% land =0
% A_1 = (land * eye(6) - A)
% solve(A_1*V ==0 ,V)
[~,~]=jordan(A);
%inv(V)*A*V
%%
s = tf('s')
syms a
C = [0 0 1 0 0 0;0 0 0 0 1 0]
D = [0;0]
G = ss2tf(A,B,C,D)
G = C*inv(s*eye(6)-A)*B+D
solve(vpa(-5*a^2 + 1.51e-13*a + 980))
solve(vpa(a^4 - 5.684e-14*a^3 - 401.8*a^2 + 1.182e-11*a + 2.113e04))
pzmap(G) % has no zeros
%%
%s = tf('s')
syms s
eig(A)
sys = ss(A,B,C,D)
w = isstable(sys) %% is not stable
fi = (inv(s*eye(6)-A))
det = det(fi)
f = fi/det
fi_t =ilaplace(f)
%%
syms t
u = 1/s
x_0 = [1;0;0;0;0;0]
X = fi*x_0 + fi*B*u
Y = C*fi*x_0 + C*fi*B*u +D *u
%%
rref(A)
[V,J]=jordan(A)
A_3 =jordan(A)
B_2 = inv(V)*B
C_2 = C*V
D_2 = D
%% 9
R=ctrb(A,B)
rank(R)
O=obsv(A,C)
rank(O)
T_2 = [0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1;1 1 0 0 0 0;0 1 0 0 0 0]
A_c = T_2*A*inv(T_2)
B_c = T_2*B
C_c = C*inv(T_2)
A_bar_o = [0 1 0 0;205.8 0 -98 0;0 0 0 1;-196 0 196 0]
A_bar_o_bar = [0 1;0 0]
A_bar_21 = [-0.98 0 0 0;-0.98 0 0 0]
A_bar_0 = [0 0;0 0;0 0;0 0]
B_bar_o = [0 -0.5 0 0]'
B_bar_o_bar = [0.5 0.5]'
C_bar_o = [1 0 0 0;0 0 1 0]
Res_3 = rank(ctrb(A_bar_o,B_bar_o))
Res_4 = rank(obsv(A_bar_o,C_bar_o))

%% 10
sys = ss(A,B,C,D);
sysr = minreal(sys);











