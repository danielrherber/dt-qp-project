%--------------------------------------------------------------------------
% DTQP_TEST_multiphase.m
% Test multiphase LQDO problem creation
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

tests = 1:2;
% tests = 1;

% go through the tests
for k = 1:length(tests)

    % test setup
    switch tests(k)
        %------------------------------------------------------------------
        case 1
        % different options for each phase
        opts.dt(1).defects = 'TR';
        opts.dt(1).nt = 105;
        opts.dt(2).mesh  = 'LGL';
        opts.dt(2).nt = 100;
        opts.dt(3).defects  = 'PS';
        opts.dt(3).quadrature  = 'G';
        opts.dt(3).mesh  = 'LGL';
        linkageflag = true;
        [setup,opts] = problem(opts,linkageflag);
        %------------------------------------------------------------------
        case 2
        % different options for each phase
        opts.dt(1).defects = 'TR';
        opts.dt(1).nt = 105;
        opts.dt(2).mesh  = 'LGL';
        opts.dt(2).nt = 100;
        opts.dt(3).defects  = 'PS';
        opts.dt(3).quadrature  = 'G';
        opts.dt(3).mesh  = 'LGL';
        linkageflag = false;
        [setup,opts] = problem(opts,linkageflag);
    end

    % run the test and time
    [T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

    % test analysis
    figure(1)
    plot(T,U);

    figure(2)
    plot(T,Y);

end

% problem structure
function [setup,opts] = problem(opts,linkageflag)

% options
opts.general.displevel = 2; % verbose

%% phase 1
phs = 1;

% time horizon
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

%% linkage equality constraints
if linkageflag
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
end

%% phase 2
phs = 2;

% time horizon
setup(phs).t0 = 1;
setup(phs).tf = 2;

% dynamics
setup(phs).A = 2*eye(2);
setup(phs).B = eye(2);

% objective
setup(phs).L.right = 1;
setup(phs).L.left = 1;
setup(phs).L.matrix = eye(2);

%% phase 3
phs = 3;

% time horizon
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

end