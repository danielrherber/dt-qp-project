%--------------------------------------------------------------------------
% SecondOrderSingular_solution.m
% Solution equations for the Second Order Singular
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

% initialize
syms t ts real % ts is the switching time
syms y1(t) y2(t) y3(t)
syms u(t)
syms y1ts y2ts y3ts real

% state derivatives
Dy1 = diff(y1,t);
Dy2 = diff(y2,t);
Dy3 = diff(y3,t);

% final time
tf = 5;

%--- symbolically determine 1st phase
% solve for states between 0 and ts with u = -1
solp1 = dsolve(Dy1==y2,Dy2==-1,Dy3==(y1^2+y2^2)/2,...
    y1(0) == 0,y2(0) == 1,y3(0) == 0);

% combine
Yp1_syms = [solp1.y1,solp1.y2,solp1.y3];

%--- symbolically determine 2nd phase
% solve for states 1 and 2 between ts and tf with u = y1
solp2_12 = dsolve(Dy1==y2,Dy2==y1,...
    y1(ts) == y1ts,y2(ts) == y2ts);

% determine Dy3 from y1 and y2
Dy3_sol2 = simplify((solp2_12.y1)^2+(solp2_12.y2)^2,'steps',100);

% solve for state 3 between ts and tf
solp2_3 = dsolve(Dy3==Dy3_sol2/2,...
    y3(ts) == y3ts);

% combine
Yp2_syms = [solp2_12.y1,solp2_12.y2,solp2_3];

%--- objective as a function of ts
% determine 1st phase
Yp1 = Yp1_syms;
Yp1ts = subs(Yp1,'t',ts);

% determine 2nd phase using 1st phase final values
obj = Yp2_syms(3);
obj = subs(obj,'y1ts',Yp1ts(1));
obj = subs(obj,'y2ts',Yp1ts(2));
obj = subs(obj,'y3ts',Yp1ts(3));
obj = subs(obj,'t',tf);

% simplify
obj = simplify(obj,'steps',100);

% create matlab function
obj = matlabFunction(obj);

%--- compute transition time
% options
options = optimset('Display','iter','TolX',eps);

% solve
[ts,F] = fminbnd(@(ts) obj(ts),1,2,options);

%--- solution with computed ts
% objective
F = sym(F);

% states - 1st phase
Yp1 = Yp1_syms;
Yp1ts = subs(Yp1,'t',ts);

% states - 2nd phase
Yp2 = Yp2_syms;
Yp2 = subs(Yp2,'y1ts',Yp1ts(1));
Yp2 = subs(Yp2,'y2ts',Yp1ts(2));
Yp2 = subs(Yp2,'y3ts',Yp1ts(3));
Yp2 = subs(Yp2,'ts',ts);

%% output the functions
path = msavename(mfilename('fullpath'),'');

prob = 'SecondOrderSingular';

aF = matlabFunction(F,'file',[path,prob,'_F']);
aU = matlabFunction(-1*heaviside(ts-t) + Yp2(1)*heaviside(t-ts),'Vars',{'t'},'file',[path,prob,'_U']);
aY = matlabFunction(Yp1*heaviside(ts-t) + Yp2*heaviside(t-ts),'Vars',{'t'},'file',[path,prob,'_Y']);

%% plot solution
% T = linspace(0,tf,1e5)';
% disp(aF)
% hf = figure; hf.Color = 'w'; hold on
% plot(T,aY(T))
% plot(T,aU(T))