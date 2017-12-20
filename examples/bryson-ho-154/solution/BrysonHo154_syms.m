clear
clc

syms t0 tf c1 c2 x0 v0 t real
syms x(t)

Dx = diff(x,t);
D2x = diff(x,t,2);

Df = (1/c2 + (tf-t)^3/3)*(1/c1 + tf -t) - (tf-t)^4/4;
Gv = (1/c2 + (tf-t)^2/c1 + (tf-t)^3/3)/Df;
Gy = ((tf-t)/c1 + (tf-t)^2/2)/Df;
u = -Gv*Dx - Gy*x;

X = dsolve( D2x == u, Dx(t0) == v0, x(t0) == x0 );

V = diff(X,t);

Xf = subs(X,t,tf);
Vf = subs(V,t,tf);

U = subs(u,x,X);
U = subs(U,Dx,V);

L = 1/2*U^2;
I = c1/2*Vf^2 + c2/2*Xf^2 + int(L,t,t0,tf,'IgnoreAnalyticConstraints',true) ;

u = simplify(U,'steps',100,'IgnoreAnalyticConstraints',true);
v = simplify(V,'steps',100,'IgnoreAnalyticConstraints',true);
x = simplify(X,'steps',100,'IgnoreAnalyticConstraints',true);
I = simplify(I,'steps',100,'IgnoreAnalyticConstraints',true);

%%
path = msavename(mfilename('fullpath'),'');

prob = 'BrysonHo154';

aU = matlabFunction(u,'file',[path,prob,'_U']);
aX = matlabFunction([x,v],'file',[path,prob,'_Y']);
aPSI = matlabFunction(I,'file',[path,prob,'_F']);