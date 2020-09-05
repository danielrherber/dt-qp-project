%--------------------------------------------------------------------------
% SimpleSASA_solution.m
% Closed-form equations for the simple SASA problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

% initialize
syms t ts % ts is the switching time
syms P1 k f J umax thetaf
syms theta(t)

Dtheta = diff(theta,t);
DDtheta = diff(Dtheta,t);

% solve for scaled theta between 0 and ts
U = 1;
theta1 = dsolve(DDtheta == -P1*theta + U, theta(0) == 0, Dtheta(0) == 0);
theta1 = simplify(theta1);

Dtheta1 = diff(theta1,t);
Dtheta1 = simplify(Dtheta1);

% solve for scaled theta between ts and tf=2*pi
U = -1;
theta2 = dsolve(DDtheta == -P1*theta + U, theta(2*pi) == thetaf, Dtheta(2*pi) == 0);
theta2 = simplify(theta2);

Dtheta2 = diff(theta2,t);
Dtheta2 = simplify(Dtheta2);

% substitute t for ts
theta1ts = subs(theta1,'t',ts);
theta2ts = subs(theta2,'t',ts);
Dtheta1ts = subs(Dtheta1,'t',ts);
Dtheta2ts = subs(Dtheta2,'t',ts);

% solve for theta at time tf
thetaf = solve( theta1ts == theta2ts, thetaf );
thetaf = simplify(thetaf,'steps',10);

% calculate gradient of the objective function wrt P1
dJdP1 = diff(thetaf,P1);
dJdP1 = simplify(dJdP1,'steps',10);
eqn1 = dJdP1 == 0; % optimality condition

% substituted thetaf for the result above to determine the continuity condition
Dtheta1ts = subs(Dtheta1ts,'thetaf',thetaf);
Dtheta2ts = subs(Dtheta2ts,'thetaf',thetaf);
eqn2 = Dtheta1ts == Dtheta2ts; % continuity condition

% numerically find the solution to the set of equations
S = vpasolve([eqn1, eqn2],[P1,ts],[0.08656, 4.578]);

% calculate scaled thetaf
thetaf = subs(thetaf,'P1',S.P1);
thetaf = subs(thetaf,'ts',S.ts);
thetaf = real(vpa(thetaf)); % sometimes a very small complex part from vpasolve

% calculate scaled theta
theta1 = subs(theta1,'P1',S.P1);
theta1 = subs(theta1,'ts',S.ts);
theta1 = subs(theta1,'thetaf',thetaf);
theta2 = subs(theta2,'P1',S.P1);
theta2 = subs(theta2,'ts',S.ts);
theta2 = subs(theta2,'thetaf',thetaf);

% create matlab functions for scaled solution
out.scaled.P1 = double(S.P1);
out.scaled.ts = double(S.ts);
out.scaled.thetaf = double(thetaf);
out.scaled.theta = matlabFunction(theta2*heaviside(t-S.ts) + theta1*heaviside(-t+S.ts),...
    'Vars',{'t'});
out.scaled.U = matlabFunction(1-2*heaviside(t-S.ts),'Vars',{'t'});

% scaling constants
ctheta = umax*f^2/(4*pi^2*J);
ct = f/(2*pi);
cu = umax;

% time substitution
theta1 = subs(theta1,t,t/ct);
theta2 = subs(theta2,t,t/ct);

% create matlab functions for unscaled solution
out.unscaled.k = matlabFunction((S.P1*4*pi^2*J)/(f^2),'Vars',{'f','J'});
out.unscaled.ts = matlabFunction(ct*S.ts,'Vars',{'f'});
out.unscaled.thetaf = matlabFunction(ctheta*thetaf,'Vars',{'umax','f','J'});
out.unscaled.theta = matlabFunction(ctheta*(theta2*heaviside(t-(ct*S.ts)) + ...
    theta1*heaviside(-t+(ct*S.ts)) ),...
    'Vars',{'t','umax','f','J'});
out.unscaled.U = matlabFunction(cu*(1-2*heaviside(t-(ct*S.ts))),'Vars',{'t','umax','f'});

%% output the functions
path = msavename(mfilename('fullpath'),'');

prob = 'SimpleSASA';
aP = matlabFunction((S.P1*4*pi^2*J)/(f^2),'Vars',{'f','J'},'file',[path,prob,'_P']);
aU = matlabFunction(cu*(1-2*heaviside(t-(ct*S.ts))),'Vars',{'t','umax','f'},'file',[path,prob,'_U']);
aF = matlabFunction(ctheta*thetaf,'Vars',{'umax','f','J'},'file',[path,prob,'_F']);
aY1 = matlabFunction(real(ctheta*(theta2*heaviside(t-(ct*S.ts)) + theta1*heaviside(-t+(ct*S.ts)))),...
    'Vars',{'t','umax','f','J'},'file',[path,prob,'_Y1']);
aY2 = matlabFunction(real(diff(ctheta*(theta2*heaviside(t-(ct*S.ts)) + theta1*heaviside(-t+(ct*S.ts))),t,1)),...
    'Vars',{'t','umax','f','J'},'file',[path,prob,'_Y2']);
aTs = matlabFunction(ct*S.ts,'Vars',{'f'},'file',[path,prob,'_Ts']);