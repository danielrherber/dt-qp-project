clear
clc

syms t real
syms t0 tf x0 real
syms a b q r m real
syms X(t) S(t)

DS = diff(S,t);
s = dsolve( DS == b^2/r*S^2 - 2*a*S - q, S(tf) == m);

s = simplify(s,'steps',1000,'IgnoreAnalyticConstraints',true);

syms c1 c2

s = subs(s,sqrt(a^2*r + b^2*q),c1);
s = subs(s,-c1*tf + r^(1/2)*atanh((- m*b^2 + a*r)/(c1*r^(1/2))),c2);

DX = diff(X,t);

x = dsolve( DX == (a - b^2/r*s)*X, X(t0) == x0);

x = simplify(x,'steps',1000,'IgnoreAnalyticConstraints',true);

syms c3
x = subs(x,x0/cosh((c2 + c1*t0)/r^(1/2)),c3);

u = -b/r*s*x;

u = simplify(u,'steps',1000,'IgnoreAnalyticConstraints',true);

L = q*x^2 + r*u^2;
xf = subs(x,'t',tf);
I = int(L,t,t0,tf) + m*xf^2;

I = simplify(I,'steps',1000,'IgnoreAnalyticConstraints',true);

C3 = x0/cosh((c2 + c1*t0)/r^(1/2));
C2 = -c1*tf + r^(1/2)*atanh((- m*b^2 + a*r)/(c1*r^(1/2)));
C1 = sqrt(a^2*r + b^2*q);

%%
path = msavename(mfilename('fullpath'),'');

prob = 'LQRScalar';

aU = matlabFunction(u,'file',[path,prob,'_U']);
aY = matlabFunction(x,'file',[path,prob,'_Y']);
aF = matlabFunction(I,'file',[path,prob,'_F']);

aC1 = matlabFunction(C1,'file',[path,prob,'_C1']);
aC2 = matlabFunction(C2,'file',[path,prob,'_C2']);
aC3 = matlabFunction(C3,'file',[path,prob,'_C3']);