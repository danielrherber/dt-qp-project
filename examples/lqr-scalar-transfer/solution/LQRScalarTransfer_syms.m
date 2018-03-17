clear; clc
syms t r q a b c d f k real
assume(f > 0)
syms Y(t) P(t)

DY = diff(Y,t);
DP = diff(P,t);

Sol = dsolve(DY==a*Y-b^2*P/(2*r),DP==-2*q*Y-a*P,Y(0)==c,Y(f)==d);
Y = Sol.Y;
P = Sol.P;
U = -b*Sol.P/(2*r);

Y = simplify(Y,100,'criterion','preferReal');
U = simplify(U,100,'criterion','preferReal');

F = int(r*U^2+q*Y^2,t,0,f);
F = simplify(F,100,'criterion','preferReal');

%% save the functions
path = msavename(mfilename('fullpath'),'');
prob = 'LQRScalarTransfer';

aU = matlabFunction(vpa(U),'file',[path,prob,'_U']);
aY = matlabFunction(Y,'file',[path,prob,'_Y']);
aF = matlabFunction(F,'file',[path,prob,'_F']);