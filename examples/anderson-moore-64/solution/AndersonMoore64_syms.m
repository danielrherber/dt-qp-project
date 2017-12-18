clear
clc

syms tf t0 x0 t real
syms x(t)

u = -0.5*(1-exp(t)*exp(-tf))/(1+exp(t)*exp(-tf))*x(t);

X = dsolve(diff(x(t),t) == 1/2*x + u, x(t0) == x0);

u = subs(u,x,X);

x = X;

L = 2*exp(-t)*u^2 + 1/2*exp(-t)*x^2;
I = int(L,t,t0,tf);

%%
path = msavename(mfilename('fullpath'),'');

prob = 'AndersonMoore64';

aU = matlabFunction(u,'file',[path,prob,'_U']);
aX = matlabFunction(x,'file',[path,prob,'_Y']);
aPSI = matlabFunction(I,'file',[path,prob,'_F']);