%--------------------------------------------------------------------------
% DTQPtest_scaling.m
% Testing scaling function
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary Contributor: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

p.ell = 1/9;

% options
opts.dt.nt = 10000;
opts.general.displevel = 0;

% time horizon 
setup.t0 = 0; setup.tf = 1;

% system dynamics
A = [0 1; 0 0]; 
B = [0;1];

% Lagrange term
L(1).left = 1; % control variables
L(1).right = 1; % control variables
L(1).matrix(1,1) = 1/2; % 1/2*u.^2

UB(1).right = 4; UB(1).matrix = [0;1]; % initial states
LB(1).right = 4; LB(1).matrix = [0;1];
UB(2).right = 5; UB(2).matrix = [0;-1]; % final states
LB(2).right = 5; LB(2).matrix = [0;-1];
UB(3).right = 2; UB(3).matrix = [p.ell;Inf]; % states

% combine structures
setup.A = A; setup.B = B; setup.L = L; setup.UB = UB; setup.LB = LB; setup.p = p;

% normal
tic
disp('normal')
[T0,U0,Y0,P0,F0] = DTQP_solve(setup,opts);
toc
disp(F0)

% scale 1
tic
disp('scale 1')
setup.scaling(1).right = 1; % controls
setup.scaling(1).matrix = [6];
[T1,U1,Y1,P1,F1] = DTQP_solve(setup,opts);
toc
disp(F1)

% scale 2
tic
disp('scale 2')
setup.scaling(1).right = 2; % states
setup.scaling(1).matrix = [1/9 1];
[T2,U2,Y2,P2,F2] = DTQP_solve(setup,opts);
toc
disp(F2)

% scale 3
tic
disp('scale 3')
setup.scaling(1).right = 2; % states
setup.scaling(1).matrix = [1/9 1];
setup.scaling(2).right = 1; % controls
setup.scaling(2).matrix = [6];
[T3,U3,Y3,P3,F3] = DTQP_solve(setup,opts);
toc
disp(F3)