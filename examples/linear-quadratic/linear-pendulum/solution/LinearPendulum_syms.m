%--------------------------------------------------------------------------
% LinearPendulum.m
% Use symoblic operations to help determine the solution for this example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
clear; close all

% initialize
syms t positive
syms x(t)
syms tf positive
syms m k umax t0 positive
syms s x0 v0 real
Dx = diff(x,t);
DDx = diff(Dx,t);

% state 1 time solution
solx = dsolve(DDx == -x + sign(s), x(t0) == x0, Dx(t0) == v0);
solx = simplify(solx,'steps',1000);

% state 2 time solution
solDx = diff(solx,t);
solDx = simplify(solDx,'steps',1000);

% control time solution
solu = umax*sign(-sin(sqrt(k/m)*(t-tf)));

% save the functions
path = msavename(mfilename('fullpath'),'');
prob = 'LinearPendulum';
ay1 = matlabFunction(solx,'file',[path,prob,'_Y1'],'Vars',[s,t0,x0,v0,t]);
ay2 = matlabFunction(solDx,'file',[path,prob,'_Y2'],'Vars',[s,t0,x0,v0,t]);
au = matlabFunction(solu,'file',[path,prob,'_U'],'Vars',[m,k,umax,tf,t]);