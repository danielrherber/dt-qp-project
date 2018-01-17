clear
clc

syms t
syms p(t) k(t) x(t)
syms b m r x0 tf real
syms a1 w1 a2 w2 real

Dx = diff(x,t);
Dp = diff(p,t);
Dk = diff(k,t);

d = a1*sin(w1*t) + a2*sin(w2*t);

p = dsolve(Dp == p^2*b^2/r, p(tf) == m);
p = simplify(p,'steps',1000,'IgnoreAnalyticConstraints',true);

k = dsolve(Dk == p*b^2/r*k - p*d, k(tf) == 0);
k = simplify(k,'steps',1000,'IgnoreAnalyticConstraints',true);

u = -(p*x + k)/r*b;

x = dsolve(Dx == b*u + d, x(0) == x0);
x = simplify(x,'steps',1000,'IgnoreAnalyticConstraints',true);

u = -(p*x + k)/r*b;
u = simplify(u,'steps',1000,'IgnoreAnalyticConstraints',true);

xtf = subs(x,'t',tf);
L = u^2;
I = int(L,t,0,tf) + m*xtf^2;
I = simplify(I,'steps',1000,'IgnoreAnalyticConstraints',true);

%%
path = msavename(mfilename('fullpath'),'');
prob = 'DTQP3';

matlabFunction(u,'file',[path,prob,'_U']);
matlabFunction(x,'file',[path,prob,'_Y']);
matlabFunction(I,'file',[path,prob,'_F']);