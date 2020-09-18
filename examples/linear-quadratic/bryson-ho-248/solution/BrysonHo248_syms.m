close all; clear; clc

syms a b c real
syms t t1 t2 tf positive
syms y1(t) y2(t)

Dy1 = diff(y1,t);
Dy2 = diff(y2,t);

% 0 < t < t1
ut1 = c;
sol1 = dsolve(Dy1 == y2 + ut1, Dy2 == -ut1, y1(0) == a, y2(0) == b);
y1t1 = simplify(subs(sol1.y1,'t',t1),'steps',100);
y2t1 = simplify(subs(sol1.y2,'t',t1),'steps',100);

% t1 < t < t2
ut2 = -(y1+y2);
sol2 = dsolve(Dy1 == y2 + ut2, Dy2 == -ut2, y1(t1) == y1t1, y2(t1) == y2t1);
y1t2 = simplify(subs(sol2.y1,'t',t2),'steps',100);
y2t2 = simplify(subs(sol2.y2,'t',t2),'steps',100);

% t2 < t < tf
ut3 = -c;
sol3 = dsolve(Dy1 == y2 + ut3, Dy2 == -ut3, y1(t2) == y1t2, y2(t2) == y2t2);
y1tf = simplify(subs(sol3.y1,'t',tf),'steps',100);
y2tf = simplify(subs(sol3.y2,'t',tf),'steps',100);

% piecewise solutions
y1 = sol3.y1*(heaviside(t-t2) - heaviside(t-tf)) + ...
    sol2.y1*(heaviside(t-t1) - heaviside(t-t2)) + ...
    sol1.y1*(heaviside(t+1) - heaviside(t-t1));

y2 = sol3.y2*(heaviside(t-t2) - heaviside(t-tf)) + ...
    sol2.y2*(heaviside(t-t1) - heaviside(t-t2)) + ...
    sol1.y2*(heaviside(t+1) - heaviside(t-t1));

u = ut3*(heaviside(t-t2) - heaviside(t-tf)) + ...
    -(sol2.y1+sol2.y2)*(heaviside(t-t1) - heaviside(t-t2)) + ...
    ut1*(heaviside(t+1) - heaviside(t-t1));

% objective function
F = 1/2*(int(sol1.y1^2,t,0,t1)+int(sol2.y1^2,t,t1,t2)+int(sol3.y1^2,t,t2,tf)); 

% implicit equations
eqs = [y1tf==0,y2tf==0];

%% save the functions
prob = 'BrysonHo248';

mypath = msavename(mfilename('fullpath'),'');

aU = matlabFunction(u,'file',[mypath,prob,'_U']);
aY = matlabFunction([y1,y2],'file',[mypath,prob,'_Y']);
aF = matlabFunction(F,'file',[mypath,prob,'_F']);

%% plot an example solution
a = -5; b = -2; c = 40; tf = 1;
t = linspace(0,tf,10000)';

[t1,t2] = BrysonHo248_T(a,b,c,tf);
U = BrysonHo248_U(a,b,c,t,t1,t2,tf);
Y = BrysonHo248_Y(a,b,c,t,t1,t2,tf);
F = BrysonHo248_F(a,b,c,t1,t2,tf);

plot(t,U); hold on
plot(t,Y); hold on

disp(F)