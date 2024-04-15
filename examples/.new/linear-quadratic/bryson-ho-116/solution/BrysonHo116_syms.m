clear
clc

syms t x0 v0 real
syms t1 t2 tf positive
syms x(t)

Dx = diff(x,t);
DDx = diff(Dx,t);

sympref('HeavisideAtOrigin', 0);
u = 1*(heaviside(t-t2) - heaviside(t-tf)) + -1*(heaviside(t+1) - heaviside(t-t1));

% 0 < t < t1
x1 = dsolve(DDx == -1, x(0) == x0, Dx(0) == v0);
Dx1 = diff(x1,t);
% t1 < t < t2
xt1 = simplify(subs(x1,'t',t1));
Dxt1 = simplify(subs(Dx1,'t',t1));
x2 = dsolve(DDx == 0, x(t1) == xt1, Dx(t1) == Dxt1);
Dx2 = diff(x2,t);
% t2 < t < tf
xt2 = simplify(subs(x2,'t',t2));
Dxt2 = simplify(subs(Dx2,'t',t2));
x3 = dsolve(DDx == 1, x(t2) == xt2, Dx(t2) == Dxt2);
Dx3 = diff(x3,t);


x = x3*(heaviside(t-t2) - heaviside(t-tf)) + ...
    x2*(heaviside(t-t1) - heaviside(t-t2)) + ...
    x1*(heaviside(t+1) - heaviside(t-t1));

x = subs(x,'t1',1/2*(tf + v0 - sqrt((tf-v0)^2-(4*x0+2*v0^2))));
x = subs(x,'t2',1/2*(tf + v0 + sqrt((tf-v0)^2-(4*x0+2*v0^2))));

Dx = Dx3*(heaviside(t-t2) - heaviside(t-tf)) + ...
    Dx2*(heaviside(t-t1) - heaviside(t-t2)) + ...
    Dx1*(heaviside(t+1) - heaviside(t-t1));

Dx = subs(Dx,'t1',1/2*(tf + v0 - sqrt((tf-v0)^2-(4*x0+2*v0^2))));
Dx = subs(Dx,'t2',1/2*(tf + v0 + sqrt((tf-v0)^2-(4*x0+2*v0^2))));

u = subs(u,'t1',1/2*(tf + v0 - sqrt((tf-v0)^2-(4*x0+2*v0^2))));
u = subs(u,'t2',1/2*(tf + v0 + sqrt((tf-v0)^2-(4*x0+2*v0^2))));

t1 = 1/2*(tf + v0 - sqrt((tf-v0)^2-(4*x0+2*v0^2)));
t2 = 1/2*(tf + v0 + sqrt((tf-v0)^2-(4*x0+2*v0^2)));

I = (t1-0) + (tf-t2);

% u = simplify(u,'steps',1000);
% x = simplify(x,'steps',1000);
% Dx = simplify(Dx,'steps',1000);
% I = simplify(I,'steps',1000);

%%
path = msavename(mfilename('fullpath'),'');

prob = 'BrysonHo116';

aU = matlabFunction(u,'file',[path,prob,'_U']);
aX = matlabFunction([x,Dx],'file',[path,prob,'_Y']);
aPSI = matlabFunction(I,'file',[path,prob,'_F']);