%--------------------------------------------------------------------------
% DTQP_TEST_multipleIntervalPS.m
% Test multiple-interval pseudospectral method approach on the
% BrysonDenham problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

tests = 1:5;
% tests = 5;

% go through the tests
for k = 1:length(tests)

    clear opts

    % test setup
    switch tests(k)
        %------------------------------------------------------------------
        case 1
        t0 = 0; tf = 1;
        ni = 20;
        xT = linspace(t0,tf,ni+1);
        opts.dt.nt = 4;
        %------------------------------------------------------------------
        case 2
        t0 = 0; tf = 1;
        ni = 4;
        xT = linspace(t0,tf,ni+1);
        opts.dt.nt = 7;
        %------------------------------------------------------------------
        case 3
        t0 = 0; tf = 1;
        ni = 3;
        xT = linspace(t0,tf,ni+1);
        opts.dt.nt = 4;
        %------------------------------------------------------------------
        case 4
        t0 = 0; tf = 1;
        ni = 100;
        xT = linspace(t0,tf,ni+1);
        opts.dt.nt = 8;
        %------------------------------------------------------------------
        case 5
        t0 = 0; tf = 1;
        ni = 4;
        xT = linspace(t0,tf,ni+1);
        opts.dt(1).nt = 4;
        opts.dt(2).nt = 5;
        opts.dt(3).nt = 6;
        opts.dt(4).nt = 10;
    end

    % run the test and time
    [setup,opts] = problem(opts,xT);
    [T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

    % test analysis
    figure(1)
    plot(T,U); hold on

    figure(2)
    plot(T,Y); hold on
end

% problem structure
function [setup,opts] = problem(opts,xT)

% tunable parameters
p.ell = 1/9;

% options
[opts.dt.defects] = deal('PS');
[opts.dt.quadrature] = deal('G');
[opts.dt.mesh] = deal('LGL');

%% setup
% system dynamics
A = [0 1;0 0]; B = [0;1];

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = 1/2; % 1/2*u^2

% simple bounds
UB(1).right = 2; UB(1).matrix = [p.ell;Inf]; % states

%% linkage equality constraints
idx = 0;
q = eye(length(A));
for k = 1:length(A)
    idx = idx + 1;
    LY(idx).left.linear.right = 5; % final states
    LY(idx).left.linear.matrix = q(:,k);
    LY(idx).right.linear.right = 4; % initial states
    LY(idx).right.linear.matrix = -q(:,k);
    LY(idx).b = 0;
end

%% phases
% combine structures
for phs = 1:length(xT)-1
    setup(phs).A = A;
    setup(phs).B = B;
    setup(phs).L = L;

    UBt = UB;
    LBt = [];

    if phs == 1
        UBt(2).right = 4; UBt(2).matrix = [0;1]; % initial states
        LBt(1).right = 4; LBt(1).matrix = [0;1];
    end

    if phs == length(xT)-1
        UBt(2).right = 5; UBt(2).matrix = [0;-1]; % final states
        LBt(1).right = 5; LBt(1).matrix = [0;-1];
    end

    setup(phs).UB = UBt;
    if ~isempty(LBt)
        setup(phs).LB = LBt;
    end

    setup(phs).t0 = xT(phs);
    setup(phs).tf = xT(phs+1);
    setup(phs).p = p;

    if phs < length(xT)-1
        setup(phs).LY = LY;
    end
end
end