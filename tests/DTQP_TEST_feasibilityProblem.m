%--------------------------------------------------------------------------
% DTQP_TEST_feasibilityProblem.m
% Solve a problem with no objective function
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

tests = 1:3;
% tests = 1;

% problem structure
setup = problem;

% options
opts.dt.nt = 100;

% go through the tests
for k = 1:length(tests)

    % test setup
    switch tests(k)
        %------------------------------------------------------------------
        case 1
        setup.tf = 6; % feasible
        %------------------------------------------------------------------
        case 2
        setup.tf = 4.17; % feasible
        %------------------------------------------------------------------
        case 3
        setup.tf = 2; % infeasible
    end

    % run the test and time
    [T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

    % test analysis
    disp(strcat("tf = ",string(setup.tf)))
    if isnan(F)
        disp('infeasible')
    else
        disp('feasible')
    end

    % plot results
    hf = figure(1); hf.Color = 'w'; hold on
    plot(T,U,'r'); plot(T,Y,'k');

end

% problem structure
function setup = problem

% parameters
g = 1.5; umax = 3;
v0 = -2; h0 = 10;

% setup
setup.t0 = 0;
setup.auxdata = [];

% system dynamics
A = [0,1;0,0];
B = [0;1];
d = [0;-g];

% linear boundary constraints
Y(1).linear(1).right = 4; % initial states
Y(1).linear(1).matrix = [1;0];
Y(1).b = h0;
Y(2).linear(1).right = 4; % initial states
Y(2).linear(1).matrix = [0;1];
Y(2).b = v0;
Y(3).linear(1).right = 5; % final states
Y(3).linear(1).matrix = [1;0];
Y(3).b = 0;
Y(4).linear(1).right = 5; % final states
Y(4).linear(1).matrix = [0;1];
Y(4).b = 0;

% simple control lower bounds
LB(1).right = 1; % control
LB(1).matrix = 0;
LB(2).right = 2; % states
LB(2).matrix = [0;-Inf];

% simple upper bounds
UB(1).right = 1; % control
UB(1).matrix = umax;

% combine
setup.A = A; setup.B = B; setup.d = d;
setup.Y = Y; setup.LB = LB; setup.UB = UB;

end