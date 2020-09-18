clear
clc

syms t real

u = -(tanh(1-t) + 0.5)*cosh(1-t)/cosh(1);
x = cosh(1-t)/cosh(1);

u = simplify(u,'steps',1000,'IgnoreAnalyticConstraints',true);
x = simplify(x,'steps',1000,'IgnoreAnalyticConstraints',true);

L = u^2 + x*u + 5/4*x^2;
I = int(L,t,0,1);

%%
path = msavename(mfilename('fullpath'),'');

prob = 'Hager1';

aU = matlabFunction(u,'file',[path,prob,'_U']);
aX = matlabFunction(x,'file',[path,prob,'_Y']);
aPSI = matlabFunction(I,'file',[path,prob,'_F']);