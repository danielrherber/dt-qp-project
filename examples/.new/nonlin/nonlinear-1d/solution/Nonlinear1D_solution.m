%--------------------------------------------------------------------------
% Nonlinear1D_solution.m
% Closed-form equations for the Nonlinear1D problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

% initialize
syms t

% constants
A = 1 - sqrt(2)*exp(5*sqrt(2));
B = -1 + exp(10*sqrt(2));

% optimal state
y1 = (A/B)^2*exp(sqrt(2)*(10-2*t))*(1 + sqrt(2)*(B-1)/(2*B-A+1)*exp(sqrt(2)*(-5+2*t)))^2;

% optimal state derivative
Dy1 = diff(y1,t);

% optimal control
u1 = (Dy1 - 2*y1)/(2*sqrt(y1));

% integrand
L = y1 + u1^2;

% objective
F = int(L,t,0,5)/2;

% problem output options
path = msavename(mfilename('fullpath'),'');
prob = 'Nonlinear1D';

% create matlab functions for unscaled solution
aU = matlabFunction(u1,'Vars',{'t'},'file',[path,prob,'_U']);
aY = matlabFunction(y1,'Vars',{'t'},'file',[path,prob,'_Y']);
aF = matlabFunction(F,'Vars',{'t'},'file',[path,prob,'_F']);