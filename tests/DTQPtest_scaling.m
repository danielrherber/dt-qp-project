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
close all
clear
clc

% time horizon and number of time points
p.t0 = 0; p.tf = 1;
p.nt = 10000;
setup.p = p;

% system dynamics
setup.A = [0 1; 0 0]; 
setup.B = [0;1];

% Lagrange term
setup.L(1).left = 1; % control variables
setup.L(1).right = 1; % control variables
setup.L(1).matrix(1,1) = 1/2; % 1/2*u.^2

% initial state conditions
setup.Y(1).linear(1).right = 4; % initial states
setup.Y(1).linear(1).matrix = [1;0];
setup.Y(1).b = 0;
setup.Y(2).linear(1).right = 4; % initial states
setup.Y(2).linear(1).matrix = [0;1];
setup.Y(2).b = 1;

% final state conditions
setup.Y(3).linear(1).right = 5; % final states
setup.Y(3).linear(1).matrix = [1;0];
setup.Y(3).b = 0;
setup.Y(4).linear(1).right = 5; % final states
setup.Y(4).linear(1).matrix = [0;1];
setup.Y(4).b = -1;

% bounds
setup.UB(1).right = 2; % states
setup.UB(1).matrix = Inf*ones(2,1);
setup.UB(1).matrix(1,1) = 1/9;

opts.displevel = 0;

tic
disp('normal')
[T0,U0,Y0,P0,F0] = DTQP_solve(setup,opts);
toc
disp(F0)

tic
disp('scale 1')
% scaling
setup.scaling(1).right = 1; % controls
setup.scaling(1).matrix = [6];
[T1,U1,Y1,P1,F1] = DTQP_solve(setup,opts);
toc
disp(F1)

tic
disp('scale 2')
% scaling
setup.scaling(1).right = 2; % states
setup.scaling(1).matrix = [1/9 1];
[T2,U2,Y2,P2,F2] = DTQP_solve(setup,opts);
toc
disp(F2)

tic
disp('scale 3')
% scaling
setup.scaling(1).right = 2; % states
setup.scaling(1).matrix = [1/9 1];
setup.scaling(2).right = 1; % controls
setup.scaling(2).matrix = [6];
[T3,U3,Y3,P3,F3] = DTQP_solve(setup,opts);
toc
disp(F3)