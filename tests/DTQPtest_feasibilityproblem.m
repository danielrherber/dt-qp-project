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
close all
clear
clc

%% tunable parameters
p.g = 1.5; p.umax = 3;
p.v0 = -2; p.h0 = 10;

p.tf = 4.17; % feasible
% p.tf = 2; % infeasible

%% setup
p.t0 = 0;
setup.p = p;

% system dynamics
setup.A = [0 1; 0 0]; 
setup.B = [0;1]; 
setup.d = [0;-p.g];

% linear boundary constraints
setup.Y(1).linear(1).right = 4; % initial states
setup.Y(1).linear(1).matrix = [1;0];
setup.Y(1).b = p.h0;
setup.Y(2).linear(1).right = 4; % initial states
setup.Y(2).linear(1).matrix = [0;1];
setup.Y(2).b = p.v0;
setup.Y(3).linear(1).right = 5; % final states
setup.Y(3).linear(1).matrix = [1;0];
setup.Y(3).b = 0;
setup.Y(4).linear(1).right = 5; % final states
setup.Y(4).linear(1).matrix = [0;1];
setup.Y(4).b = 0;

% linear inequality constraint 
% setup.Z(1).linear.right = 2; % states 
% setup.Z(1).linear.matrix = [-1.1 -0.4]';
% setup.Z(1).b = 0;

% simple control lower bounds
setup.LB(1).right = 1; % control
setup.LB(1).matrix = 0;
setup.LB(2).right = 2; % states
setup.LB(2).matrix = [0;-Inf];

% simple upper bounds
setup.UB(1).right = 1; % control
setup.UB(1).matrix = p.umax;

%% solve
[T,U,Y,P,F,p,opts] = DTQP_solve(setup,[]);

%% plots
if isnan(F)
    disp('infeasible')
else
    disp('feasible')
end
figure
plot(T,[Y,U])