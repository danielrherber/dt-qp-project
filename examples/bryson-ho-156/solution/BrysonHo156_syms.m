clear
clc

syms t0 tf x0 v0 w c t real
syms x(t)

Dx = diff(x,t);
D2x = diff(x,t,2);

num = (4*w^2*cos(w*(tf-t)))*x(t) + (4*w*sin(w*(tf-t)))*Dx;
den = 4*w^3/c + 2*w*(tf-t) - sin(2*w*(tf-t));
u2 = -num/den*sin(w*(tf-t));
num = (4*w^2*cos(w*(tf-t0)))*x0 + (4*w*sin(w*(tf-t0)))*v0;
den = 4*w^3/c + 2*w*(tf-t0) - sin(2*w*(tf-t0));
u = -num/den*sin(w*(tf-t));

x = dsolve(D2x == -w^2*x + u, x(t0) == x0, Dx(t0) == v0);

v = diff(x,t);

%
num = x0*cos(w*(tf-t0)) + v0*sin(w*(tf-t0))/w;
den = 1 + c/4*w^3*(2*w*(tf-t0) - sin(2*w*(tf-t0)));
xf = num/den;

L = 1/2*u^2;
I = c/2*xf^2 + int(L,t,t0,tf,'IgnoreAnalyticConstraints',true) ;

u = simplify(u,'steps',100,'IgnoreAnalyticConstraints',true);
v = simplify(v,'steps',100,'IgnoreAnalyticConstraints',true);
x = simplify(x,'steps',100,'IgnoreAnalyticConstraints',true);
I = simplify(I,'steps',100,'IgnoreAnalyticConstraints',true);

%%
path = msavename(mfilename('fullpath'),'');

prob = 'BrysonHo156';

aU = matlabFunction(u,'file',[path,prob,'_U']);
aX = matlabFunction([x,v],'file',[path,prob,'_Y']);
aPSI = matlabFunction(I,'file',[path,prob,'_F']);