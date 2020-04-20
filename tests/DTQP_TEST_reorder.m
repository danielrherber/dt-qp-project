%--------------------------------------------------------------------------
% DTQP_TEST_reorder.m
% Test reordering function
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

% problem structure
setup = problem;

% options
opts.dt.nt = 10000;
opts.general.displevel = 0;

% warmup
opts.qp.reorder = 1;
[~,~,~,~,~] = DTQP_solve(setup,opts);

% normal
tic
disp('normal')
opts.qp.reorder = 0;
[T1,U1,Y1,P1,F1,~,~] = DTQP_solve(setup,opts);
toc

% reordered
tic
disp('reordered')
opts.qp.reorder = 1;
[T2,U2,Y2,P2,F2,~,~] = DTQP_solve(setup,opts);
toc

% display different in objective function values
disp(abs(F1-F2))

% difference between state and control solutions
figure; hold on
plot(T1,U1-U2)
plot(T1,Y1-Y2)

% problem structure
function setup = problem
p.ell = 1/9;

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