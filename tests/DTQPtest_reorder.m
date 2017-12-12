%--------------------------------------------------------------------------
% DTQPtest_reorder.m
% Test reordering function
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all
clear
clc

p.t0 = 0; p.tf = 1;
p.nt = 100;
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

% warmup
opts.reorder = 1;
[~,~,~,~,~] = DTQP_solve(setup,opts);

tic
disp('normal')
opts.reorder = 0;
[T1,U1,Y1,P1,F1,~,~] = DTQP_solve(setup,opts);
toc

tic
disp('reordered')
opts.reorder = 1;
[T2,U2,Y2,P2,F2,~,~] = DTQP_solve(setup,opts);
toc

disp(abs(F1-F2))
