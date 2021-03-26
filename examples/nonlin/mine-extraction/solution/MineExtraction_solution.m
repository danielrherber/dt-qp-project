%--------------------------------------------------------------------------
% MineExtraction_solution.m
% Closed-form equations for the Mine Extraction problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

% initialize
syms t a T x0 real

% state
y = (4 - a*t + a*T)^2/(4 + a*T)^2*x0;

% control
u = -diff(y,t);

% integrand
L = -a*u + u^2/y;
L = simplify(L,'steps',10);

% objective
I = int(L,t,0,T);

%% output the functions
path = msavename(mfilename('fullpath'),'');

prob = 'MineExtraction';
aY = matlabFunction(y,'Vars',{'t','a','T','x0'},'file',[path,prob,'_Y']);
aU = matlabFunction(u,'Vars',{'t','a','T','x0'},'file',[path,prob,'_U']);
aF = matlabFunction(I,'Vars',{'t','a','T','x0'},'file',[path,prob,'_F']);