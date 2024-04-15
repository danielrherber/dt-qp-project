clear
clc

syms t0 tf x0 c t real
syms x(t)

u = -x0/(1/c + tf - t0);

xf = x0/(1 + c*(tf-t0));
x = -c*xf*(t-t0) + x0;

u = simplify(u,'steps',1000,'IgnoreAnalyticConstraints',true);
x = simplify(x,'steps',1000,'IgnoreAnalyticConstraints',true);


L = 1/2*u^2;
I = int(L,t,t0,tf) + 1/2*c*xf^2;

%%
path = msavename(mfilename('fullpath'),'');

prob = 'BrysonHo153';

aU = matlabFunction(u,'file',[path,prob,'_U']);
aX = matlabFunction(x,'file',[path,prob,'_Y']);
aPSI = matlabFunction(I,'file',[path,prob,'_F']);