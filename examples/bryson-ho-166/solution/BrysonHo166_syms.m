clear
clc

syms tf x0 v0 t real
syms x(t)

u = -2/(tf^2-(sin(tf))^2)*[sin(tf-t)*sin(tf) - tf*sin(t), -cos(tf-t)*sin(tf) + tf*cos(t)]*[x0;v0];

u = simplify(u,'steps',1000,'IgnoreAnalyticConstraints',true);

Dx = diff(x,t);
D2x = diff(x,t,2);

x = dsolve(D2x == -x + u, x(0) == x0, Dx(0) == v0);

x = simplify(x,'steps',1000,'IgnoreAnalyticConstraints',true);

v = diff(x,t);

v = simplify(v,'steps',1000,'IgnoreAnalyticConstraints',true);

L = 1/2*u^2;
I = int(L,t,0,tf);

%%
path = msavename(mfilename('fullpath'),'');

prob = 'BrysonHo166';

aU = matlabFunction(u,'file',[path,prob,'_U']);
aX = matlabFunction([x,v],'file',[path,prob,'_X']);
aPSI = matlabFunction(I,'file',[path,prob,'_PSI']);