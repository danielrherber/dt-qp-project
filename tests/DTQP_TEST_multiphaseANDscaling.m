%--------------------------------------------------------------------------
% DTQP_TEST_multiphaseANDscaling.m
% Test multiphase LQDO problem creation with scaling in each phase
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

scaleflag = true;

% options
opts.general.displevel = 2; % verbose

% different options for each phase
opts.dt(1).defects = 'TR';
opts.dt(1).nt = 105;

opts.dt(2).mesh  = 'LGL';
opts.dt(2).nt = 100;

opts.dt(3).defects  = 'PS';
opts.dt(3).quadrature  = 'G';
opts.dt(3).mesh  = 'LGL';

%% phase 1
phs = 1;

setup(phs).t0 = 0;
setup(phs).tf = 1;

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

% scaling
if scaleflag
    setup(phs).scaling(1).right = 1; % controls
    setup(phs).scaling(1).matrix = [1 0.5];
    setup(phs).scaling(1).right = 2; % states
    setup(phs).scaling(1).matrix = [1 1];
end

%% linkage equality constraints
idx = 0;
q = eye(length(setup(1).A));
for k = 1:length(setup(1).A)
    idx = idx + 1;
    LY(idx).left.linear.right = 5; % final states
    LY(idx).left.linear.matrix = q(:,k);
    LY(idx).right.linear.right = 4; % initial states
    LY(idx).right.linear.matrix = -q(:,k);
    LY(idx).b = 0;
end
setup(1).LY = LY;
setup(2).LY = LY;

%% phase 2
phs = 2;

setup(phs).t0 = 1;
setup(phs).tf = 2;

% dynamics
setup(phs).A = 2*eye(2);
setup(phs).B = eye(2);

% objective
setup(phs).L.right = 1;
setup(phs).L.left = 1;
setup(phs).L.matrix = eye(2);

% scaling
if scaleflag
    setup(phs).scaling(1).right = 1; % controls
    setup(phs).scaling(1).matrix = [1 0.5];
    setup(phs).scaling(1).right = 2; % states
    setup(phs).scaling(1).matrix = [2 5];
end

%% phase 3
phs = 3;

setup(phs).t0 = 2;
setup(phs).tf = 3;

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

% scaling
if scaleflag
    setup(phs).scaling(1).right = 1; % controls
    setup(phs).scaling(1).matrix = [1 0.5];
    setup(phs).scaling(1).right = 2; % states
    setup(phs).scaling(1).matrix = [2 5];
end

%% solve
[T,U,Y,P,F,p,opts] = DTQP_solve(setup,opts);

%% plot
figure
plot(T,U);
title("controls"); legend(["control 1","control 2"])

figure
plot(T,Y);
title("states"); legend(["state 1","state 2"])