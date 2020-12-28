%--------------------------------------------------------------------------
% HagerHouRao1_solution.m
% Closed-form equations for the HagerHouRao1 problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

% initialize
syms t

% optimal controls
u1 = 2*(exp(3*t) - exp(3))/((2 + exp(3))*exp(3*t/2));
u2 = -cosh(1-t)*(tanh(1-t) + 0.5)/cosh(1);
U = [u1,u2];

% optimal states
y1 = cosh(1-t)*(2*exp(3*t) + exp(3))/((2 + exp(3))*exp(3*t/2)*cosh(1));
y2 = cosh(1)/cosh(1-t);
Y = [y1,y2];

% integrand
L = 2*y1^2*y2^2 + 1.25/(y2^2) + u2/y2 + u1^2 + u2^2;

% simplify
L = simplify(L);

% objective
F = int(L,t,0,1);

% problem output options
path = msavename(mfilename('fullpath'),'');
prob = 'HagerHouRao1';

% create matlab functions for unscaled solution
aU = matlabFunction(U,'Vars',{'t'},'file',[path,prob,'_U']);
aY = matlabFunction(Y,'Vars',{'t'},'file',[path,prob,'_Y']);
aF = matlabFunction(F,'Vars',{'t'},'file',[path,prob,'_F']);