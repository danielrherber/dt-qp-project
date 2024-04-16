clear
clc

syms tf t real
syms x(t)

u = 1 - exp(t-tf);

Dx = diff(x);
D2x = diff(x,2);
x = dsolve(D2x == -Dx + u, x(0) == 0, Dx(0) == 0);

v = diff(x,t);

u = simplify(u,'steps',1000,'IgnoreAnalyticConstraints',true);
v = simplify(v,'steps',1000,'IgnoreAnalyticConstraints',true);
x = simplify(x,'steps',1000,'IgnoreAnalyticConstraints',true);

X = matlabFunction(x);

L = 1/2*u^2;
I = int(L,t,0,tf) - X(tf,tf);

%%
path = msavename(mfilename('fullpath'),'');

prob = 'Cart';

aU = matlabFunction(u,'file',[path,prob,'_U']);
aX = matlabFunction([x,v],'file',[path,prob,'_X']);
aPSI = matlabFunction(I,'file',[path,prob,'_PSI']);