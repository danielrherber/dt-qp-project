clear
clc

syms t
syms p(t) k(t) x(t)
syms a b m q r x0 tf d0 h w real
syms h1 w real
syms A B M Q R
 
Dx = diff(x,t);
Dp = diff(p,t);
Dk = diff(k,t);

a = A;
b = B*sin(w*t);
m = M;
q = 0;
r = R;

p = dsolve(Dp == -2*p*a + p^2*b^2/r - q, p(tf) == m);
p = simplify(p,'steps',100,'IgnoreAnalyticConstraints',true);

x = dsolve(Dx == (a - p/r*b^2)*x, x(0) == x0);
x = simplify(x,'steps',100,'IgnoreAnalyticConstraints',true);

u = -(p*x)/r*b;
u = simplify(u,'steps',100,'IgnoreAnalyticConstraints',true);

xtf = subs(x,'t',tf);

L = u^2;
L = simplify(L,'steps',100,'IgnoreAnalyticConstraints',true);
% I = int(L,t,0,tf) + m*xtf^2;
% I = simplify(I,'steps',100,'IgnoreAnalyticConstraints',true);

%%
path = msavename(mfilename('fullpath'),'');
prob = 'DTQP2';

matlabFunction(u,'file',[path,prob,'_U']);
matlabFunction(x,'file',[path,prob,'_Y']);
matlabFunction(L,'file',[path,prob,'_L']);

% %%
figure
A = 1; B = 1; M = 1; R = 1; tf = 15; x0 = 1; w1 = pi;
t = linspace(0,tf,10000);
plot(t,DTQP2_U(A,B,M,R,t,tf,w1,x0)); hold on
plot(t,DTQP2_Y(A,B,M,R,t,tf,w1,x0)); hold on
