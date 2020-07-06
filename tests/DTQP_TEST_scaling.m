%--------------------------------------------------------------------------
% DTQP_TEST_scaling.m
% Testing scaling functionality
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

% tests to run (see below)
tests = [-12:12]; % all tests
% tests = 0:12; % all tests (with row scaling)
% tests = [0,1:7]; % scalar tests
% tests = [0,8:9]; % time-varying function tests
% tests = [0,10:12]; % time-based matrix tests
% tests = [-1,1]; % row scaling enabled and disabled

% problem structure
[setup,opts] = problem;

% number of tests
ntests = length(tests);

% initialize
[T,U,Y,P,F] = deal(cell(ntests,1));

% go through each test
for k = 1:length(tests)

    rng(6783097)

    % potentially remove the scaling field from previous tests
    if isfield(setup,'scaling')
        setup = rmfield(setup,'scaling');
    end

    if tests(k) > 0
        opts.method.scalematrixrows = true; % enabled
    else
        opts.method.scalematrixrows = false; % disabled
    end

    % test setup
    switch abs(tests(k))
        case 0 % no scaling
        %------------------------------------------------------------------
        case 1 % scalar (controls)
        setup.scaling(1).right = 1; % controls
        setup.scaling(1).matrix = [6];
        %------------------------------------------------------------------
        case 2 % scalar (states)
        setup.scaling(1).right = 2; % states
        setup.scaling(1).matrix = [1/9 1];
        %------------------------------------------------------------------
        case 3 % scalar (states and controls)
        setup.scaling(1).right = 2; % states
        setup.scaling(1).matrix = [1/9 1];
        setup.scaling(2).right = 1; % controls
        setup.scaling(2).matrix = [6];
        %------------------------------------------------------------------
        case 4 % scalar (controls)
        setup.scaling(1).right = 1; % controls
        setup.scaling(1).matrix = [1];
        setup.scaling(1).constant = [4];
        %------------------------------------------------------------------
        case 5 % scalar (controls)
        setup.scaling(1).right = 1; % controls
        setup.scaling(1).matrix = [pi*1e-7];
        setup.scaling(1).constant = [-exp(1)];
        %------------------------------------------------------------------
        case 6 % scalar (states)
        setup.scaling(1).right = 2; % states
        setup.scaling(1).matrix = [1/9 pi];
        setup.scaling(1).constant = [-exp(1) 7];
        %------------------------------------------------------------------
        case 7 % scalar (states and controls)
        setup.scaling(1).right = 1; % controls
        setup.scaling(1).matrix = [2];
        setup.scaling(1).constant = [-0.5];
        setup.scaling(2).right = 2; % states
        setup.scaling(2).matrix = [1/9 pi];
        setup.scaling(2).constant = [-exp(1) 7];
        %------------------------------------------------------------------
        case 8 % time-varying function (controls)
        Told = linspace(setup.t0,setup.tf,10);
        Uold = rand(length(Told),1);
        setup.scaling(1).right = 1; % controls
        setup.scaling(1).matrix = @(t) interp1(Told,Uold,t);
        setup.scaling(1).constant = @(t) interp1(Told,Uold,t);
        %------------------------------------------------------------------
        case 9 % time-varying function (controls and states)
        Told = linspace(setup.t0,setup.tf,10000);
        Uold = BrysonDenham_U(Told,setup.p.ell);
        Yold = BrysonDenham_Y(Told,setup.p.ell);
        setup.scaling(1).right = 1; % controls
        setup.scaling(1).constant = @(t) interp1(Told,Uold,t);
        setup.scaling(2).right = 2; % states
        setup.scaling(2).constant = @(t) interp1(Told,Yold,t);
        %------------------------------------------------------------------
        case 10 % time-based matrix (controls)
        setup.scaling(1).right = 1; % controls
        setup.scaling(1).matrix = rand(opts.dt.nt,1);
        setup.scaling(1).constant = -10000*rand(opts.dt.nt,1);
        %------------------------------------------------------------------
        case 11 % time-based matrix (controls and states)
        Told = linspace(setup.t0,setup.tf,opts.dt.nt);
        Uold = BrysonDenham_U(Told,setup.p.ell);
        Yold = BrysonDenham_Y(Told,setup.p.ell);
        setup.scaling(1).right = 1; % controls
        setup.scaling(1).constant = Uold;
        setup.scaling(2).right = 2; % states
        setup.scaling(2).constant = Yold;
        %------------------------------------------------------------------
        case 12 % time-based matrix (controls and states)
        Told = linspace(setup.t0,setup.tf,opts.dt.nt);
        setup.scaling(1).right = 1; % controls
        setup.scaling(1).matrix = 100*rand(opts.dt.nt,1);
        setup.scaling(1).constant = 100*rand(opts.dt.nt,1);
        setup.scaling(2).right = 2; % states
        setup.scaling(2).matrix = 100*rand(opts.dt.nt,2);
        setup.scaling(2).constant = -100*rand(opts.dt.nt,2);
    end

    % run the test and time
    t1 = tic;
    [T{k},U{k},Y{k},P{k},F{k},in(k),opts2(k)] = DTQP_solve(setup,opts);
    toc(t1)

    % test analysis
    disp(strcat("Test #",string(tests(k))))
    disp(strcat("F: ",string(F{k})))

end

%% figures
% state 1
hf = figure; hold on; hf.Color = 'w';
legendstr = strings(ntests,1);
for k = 1:ntests
    Yactual = BrysonDenham_Y(T{k},setup.p.ell);
    d = abs(Y{k}(:,1)-Yactual(:,1));
    plot(T{k},d,'linewidth',2);
    legendstr(k) = string(tests(k));
    EY1(k) = max(d);
end
legend(strcat("test ",legendstr))
xlabel("t"); ylabel("Y1 error")

% state 2
hf = figure; hold on; hf.Color = 'w';
legendstr = strings(ntests,1);
for k = 1:ntests
    Yactual = BrysonDenham_Y(T{k},setup.p.ell);
    d = abs(Y{k}(:,2)-Yactual(:,2));
    plot(T{k},d,'linewidth',2);
    legendstr(k) = string(tests(k));
    EY2(k) = max(d);
end
legend(strcat("test ",legendstr))
xlabel("t"); ylabel("Y2 error")

% control 1
hf = figure; hold on; hf.Color = 'w';
legendstr = strings(ntests,1);
for k = 1:ntests
    Uactual = BrysonDenham_U(T{k},setup.p.ell);
    d = abs(U{k}(:,1)-Uactual(:,1));
    plot(T{k},d,'linewidth',2);
    legendstr(k) = string(tests(k));
    EU(k) = max(d);
end
legend(strcat("test ",legendstr))
xlabel("t"); ylabel("U1 error")

% qp solving time vs. errors
hf = figure; hold on; hf.Color = 'w';
ha = gca; ha.YScale = 'log';
c = lines(3);
for k = 1:ntests
    if tests(k) > 0
        marker = '.';
    else
        marker = '*';
    end
    Timer = opts2(k).timer.qpsolver;
    plot(Timer,EY1(k),marker,'color',c(1,:))
    plot(Timer,EY2(k),marker,'color',c(2,:))
    plot(Timer,EU(k),marker,'color',c(3,:))
end
xlabel("timer (s)"); ylabel("error")

return

% problem structure
function [setup,opts] = problem

% BrysonDenham path constraint parameter
p.ell = 1/9;

% options
opts.dt.nt = 1000; % 10000
opts.general.displevel = 1;

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

end