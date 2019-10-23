close all; clear; clc

% initialize
syms t 
syms Y1(t) Y2(t)
syms L1(t) L2(t)
Y = [Y1;Y2];
dYdt = diff(Y,t);
L  = [L1;L2];
dLdt = diff(L,t);

% problem parameters
A = [0,1;0,-1];
B = [0;1];
Q = eye(2);
R = 0.005;
y0 = [0;-1];

% control
u = -inv(R)*B'*L;

%  state dynamics
eq1 = dYdt == A*Y + B*u;

% costate dynamics
eq2 = -dLdt == Q*Y + A'*L;

% solve
sol = dsolve(eq1,eq2,L1(1) == 0,L2(1) == 0,Y1(0) == y0(1),Y2(0) == y0(2));

% extract and simplify
Y1 = sol.Y1;
Y1 = simplify(Y1,'steps',100);
Y2 = sol.Y2;
Y2 = simplify(Y2,'steps',100);
L1 = sol.L1;
L21 = simplify(L1,'steps',100);
L2 = sol.L2;
L2 = simplify(L2,'steps',100);

% combine
Y = [Y1;Y2];
L = [L1;L2];

% optimal control
U = -inv(R)*B'*L;
U = simplify(U,'steps',100);

% integrand
I1 = simplify(expand(Y1^2),'steps',1);
I2 = simplify(expand(Y2^2),'steps',1);
I3 = simplify(expand(R*U^2),'steps',1);
I = I1 + I2 + I3;

% convert
fI = matlabFunction(I);

% numerically compute objective function value(symbolic is too expensive)
F = integral(@(t) fI(t),0,1,'AbsTol',eps,'RelTol',eps);

%% save solution
path = msavename(mfilename('fullpath'),'');
prob = 'JadduShimemura0';

% output functions
aU = matlabFunction(U,'file',[path,prob,'_U']);
aY = matlabFunction([Y1,Y2],'file',[path,prob,'_Y']);
aF = matlabFunction(sym(F),'file',[path,prob,'_F']);

%% plot the solution (for debugging)
T = linspace(0,1,10000)';

figure; hold on
plot(T,aY(T))

figure; hold on
plot(T,aU(T))