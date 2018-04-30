%--------------------------------------------------------------------------
% DTQPtest_feasibilityproblem.m
% Solve a problem with no objective function
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

opts.dt.nt = 100;

%% tunable parameters
g = 1.5; umax = 3;
v0 = -2; h0 = 10;

% p.tf = 6; % feasible
p.tf = 4.17; % feasible
% p.tf = 2; % infeasible

%% setup
p.t0 = 0;
setup.p = p;

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

% linear inequality constraint 
% setup.Z(1).linear.right = 2; % states 
% setup.Z(1).linear.matrix = [-1.1 -0.4]';
% setup.Z(1).b = 0;

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
setup.Y = Y;
setup.LB = LB;
setup.UB = UB;

%% solve
[T,U,Y,P,F,p,opts] = DTQP_solve(setup,opts);

%% plots
if isnan(F)
    disp('infeasible')
else
    disp('feasible')
end
figure
plot(T,[Y,U])