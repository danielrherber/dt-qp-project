%--------------------------------------------------------------------------
% DTQP_TEST_scaling_nonlin.m
% Testing scaling functionality for nonlinear problems
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

% tests to run (see below)
tests = [0:8];

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
        %------------------------------------------------------------------
        case 0 % no scaling
        %------------------------------------------------------------------
        case 1 % scalar (controls)
        setup.scaling(1).right = 1; % controls
        setup.scaling(1).matrix = [75];
        %------------------------------------------------------------------
        case 2 % scalar (states)
        setup.scaling(1).right = 2; % states
        setup.scaling(1).matrix = [8.6288e+03 4.3144e+03 15];
        %------------------------------------------------------------------
        case 3 % scalar (controls and states)
        setup.scaling(1).right = 1; % controls
        setup.scaling(1).matrix = [75];
        setup.scaling(2).right = 2; % states
        setup.scaling(2).matrix = [8.6288e+03 4.3144e+03 15];
        %------------------------------------------------------------------
        case 4 % scalar (controls)
        setup.scaling(1).right = 1; % controls
        setup.scaling(1).matrix = [75];
        setup.scaling(1).constant = [45];
        %------------------------------------------------------------------
        case 5 % scalar (states)
        setup.scaling(1).right = 2; % states
        setup.scaling(1).matrix = [8.6288e+03 4.3144e+03 15];
        setup.scaling(1).constant = [-10 1000 0];
        %------------------------------------------------------------------
        case 6 % scalar (controls and states)
        setup.scaling(1).right = 1; % controls
        setup.scaling(1).matrix = [75];
        setup.scaling(1).constant = [45];
        setup.scaling(2).right = 2; % states
        setup.scaling(2).matrix = [8.6288e+03 4.3144e+03 15];
        setup.scaling(2).constant = [-10 1000 0];
        %------------------------------------------------------------------
        case 7 % time-varying function (controls)
        Told = linspace(setup.t0,setup.tf,10);
        Uold = 100*rand(length(Told),1);
        setup.scaling(1).right = 1; % controls
        setup.scaling(1).matrix = @(t) interp1(Told,Uold,t);
        setup.scaling(1).constant = @(t) interp1(Told,Uold,t);
        %------------------------------------------------------------------
        case 8 % time-based matrix (controls and states)
        Told = linspace(setup.t0,setup.tf,opts.dt.nt);
        setup.scaling(1).right = 1; % controls
        setup.scaling(1).matrix = 100*rand(opts.dt.nt,1);
        setup.scaling(1).constant = 100*rand(opts.dt.nt,1);
        setup.scaling(2).right = 2; % states
        setup.scaling(2).matrix = 100*ones(opts.dt.nt,3);
        setup.scaling(2).constant = -100*ones(opts.dt.nt,3);
        %------------------------------------------------------------------
        case 9 % (NOT WORKING?) time-based matrix states
        Told = linspace(setup.t0,setup.tf,opts.dt.nt);
        setup.scaling(1).right = 2; % states
        setup.scaling(1).matrix = rand(opts.dt.nt,3);
        %------------------------------------------------------------------
    end

    % run the test and time
    t1 = tic;
    [T{k},U{k},Y{k},P{k},F{k},in2(k),opts2(k)] = DTQP_solve(setup,opts);


    % test analysis
    disp(strcat("Test #",string(tests(k))))
    disp(strcat("F: ",string(F{k})))
    toc(t1)

end

%% figures
% state 1
hf = figure; hold on; hf.Color = 'w';
legendstr = strings(ntests,1);
for k = 1:ntests
    d = abs(Y{k}(:,1));
    plot(T{k},d,'linewidth',2);
    legendstr(k) = string(tests(k));
end
ha = gca; ha.YScale = 'log';
legend(strcat("test ",legendstr))
xlabel("t"); ylabel("Y1")

% state 2
hf = figure; hold on; hf.Color = 'w';
legendstr = strings(ntests,1);
for k = 1:ntests
    d = abs(Y{k}(:,2));
    plot(T{k},d,'linewidth',2);
    legendstr(k) = string(tests(k));
end
ha = gca; ha.YScale = 'log';
legend(strcat("test ",legendstr))
xlabel("t"); ylabel("Y2")

% state 3
hf = figure; hold on; hf.Color = 'w';
legendstr = strings(ntests,1);
for k = 1:ntests
    d = abs(Y{k}(:,3));
    plot(T{k},d,'linewidth',2);
    legendstr(k) = string(tests(k));
end
ha = gca; ha.YScale = 'log';
legend(strcat("test ",legendstr))
xlabel("t"); ylabel("Y3")

% control 1
hf = figure; hold on; hf.Color = 'w';
legendstr = strings(ntests,1);
for k = 1:ntests
    d = abs(U{k}(:,1));
    plot(T{k},d,'linewidth',2);
    legendstr(k) = string(tests(k));
end
ha = gca; ha.YScale = 'log';
legend(strcat("test ",legendstr))
xlabel("t"); ylabel("U1")

return

% problem structure
function [setup,opts] = problem

%% tunable parameters
A0 = 15;
b = 5.85; % per day
Mew = 0.02; % per day
G = 0.15; % per mg of dose per day
zeta = 0.084; % per day
D = 0.00873; % per mm^2 per day
p0 = (((b-Mew)/D)^(3/2))/2;
q0 = p0/2;
umax = 75;
y0 = [p0,q0,0];
ymin = [0.1,0.1,-inf];
auxdata.y0 = y0;

%% setup
% time horizon
auxdata.t0 = 0; auxdata.tf = 1.2;

% system dynamics
element.dynamics = '[-zeta*y1*log(y1/y2); y2*(b-(Mew+(D*(y1^(2/3)))+G*u1)); u1]';
element.parameter_list = 'zeta b Mew D G';
element.parameter_values = [zeta b Mew D G];

n.ny = 3;
n.nu = 1;
setup.n = n;

% Mayer term
M(1).right = 5; % final states
M(1).left = 0; % singleton
M(1).matrix = [1,0,0];

% simple bounds
UB(1).right = 4; UB(1).matrix = y0'; % initial states
LB(1).right = 4; LB(1).matrix = y0';
UB(2).right = 1; UB(2).matrix = umax; % controls
LB(2).right = 1; LB(2).matrix = 0;
UB(3).right = 5; UB(3).matrix = [inf,inf,A0]'; % final states
LB(3).right = 2; LB(3).matrix = ymin'; % states

% guess
Y0 = [[y0];[y0]];
U0 = [[0.5*umax];[0.5*umax]];
setup.guess.X = [U0,Y0];

% combine structures
setup.element = element; setup.M = M; setup.UB = UB; setup.LB = LB;
setup.t0 = auxdata.t0; setup.tf = auxdata.tf; setup.auxdata = auxdata;

% options
opts.general.displevel = 1;
opts.general.plotflag = 1;
opts.dt.defects = 'TR';
opts.dt.quadrature = 'CTR';
opts.dt.mesh = 'ED';
opts.dt.nt = 50; % number of nodes
opts.solver.display = 'none'; % iterations
opts.solver.function = 'ipfmincon';
opts.solver.tolerance = 1e-8;
opts.method.form = 'nonlinearprogram';
opts.solver.maxiters = 4000;
opts.method.olqflag = true;

end