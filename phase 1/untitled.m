clc;
clear all;
syms m1 m2 m3 g r1 r2
syms x1 x2 x3 x4 x5 x6
a = m2+m3
b = cos(x3)*cos(x5)
c = sin(x3-x5)
d = m2*a*sin(x3)^2 + m1*m3*c^2 + m1*m2
e = 2*m2^2 + m1*m3 + 2*m2*m3 + 2*m1*m2 + m1*m3*cos(2*x5)
f = m2*a*sin(2*x3) - m1*m3*sin(2*(x5-x3))
h = m2*m3*sin(x3)*cos(x5) + m1*m3*c
Q=(m2*g*a*sin(x3)*cos(x3))/(-d*x3)
W=m2*r1*a*x4*sin(x3)/d
R=m2*m3*r2*x6*(sin(x5-2*x3)-sin(x5))/(-2*d)
T=g*e*sin(x3)/(2*r1*d*x3)
Y=f*x4/(-2*d)
U=m1*m3*g*b*sin(x5)/(-r1*d*x5)
I=r2*h*x6/(-r1*d)
O=m1*g*a*b*sin(x3)/(-r2*d*x3)
P=m1*r1*a*c*x4/(r2*d)
S=m1*g*a*((cos(x3))^2)*sin(x5)/(r2*d*x5)
D=m1*m3*x6*sin(2*(x5-x3))/(-2*d)
A_n = [0 1 0 0 0 0;0 0 Q W 0 R;0 0 0 1 0 0;0 0 T Y U I;0 0 0 0 0 1;0 0 O P S D] 