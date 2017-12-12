%--------------------------------------------------------------------------
% DTQPtest_multiphase.m
% Test multiphase LQDO problem creation
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

% options
opts.displevel = 2; % verbose

%% phase 1
phs = 1;

setup(phs).p.t0 = 0;
setup(phs).p.tf = 1;

% dynamics
setup(phs).A = -eye(2);
setup(phs).B = eye(2);

% objective
setup(phs).L.right = 1;
setup(phs).L.left = 1;
setup(phs).L.matrix = eye(2);

% initial conditions
setup(phs).LB.right = 4;
setup(phs).LB.matrix = [1;1];
setup(phs).UB.right = 4;
setup(phs).UB.matrix = [1;1];

%% phase 2
phs = 2;

setup(phs).p.t0 = 1;
setup(phs).p.tf = 2;
setup(phs).p.nt = 150;

% dynamics
setup(phs).A = 2*eye(2);
setup(phs).B = eye(2);

% objective
setup(phs).L.right = 1;
setup(phs).L.left = 1;
setup(phs).L.matrix = eye(2);

%% phase 3
phs = 3;

setup(phs).p.t0 = 2;
setup(phs).p.tf = 3;

% dynamics
setup(phs).A = -eye(2);
setup(phs).B = eye(2);

% objective
setup(phs).L.right = 1;
setup(phs).L.left = 1;
setup(phs).L.matrix = eye(2);

% final conditions
setup(phs).LB.right = 5;
setup(phs).LB.matrix = [-1;2];
setup(phs).UB.right = 5;
setup(phs).UB.matrix = [-1;2];

%% solve
[T,U,Y,P,F,p,opts] = DTQP_multiphase(setup,opts);

%%  plot
figure
plot(T,U);

figure
plot(T,Y);