%--------------------------------------------------------------------------
% DTQP_TEST_qlin_parameters.m
% A collection of test cases using qlin with parameters
% Comparisons made to SQP with fmincon
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
clc; clear; close all

testnum = 3;

% system dynamics
symb.D = '[0]';

% Lagrange term
switch testnum
    case 1
        symb.Ob = '(p1-1)^4 + p1^2 + (p2-pi)^2';
    case 2
        symb.Ob = '(p1-1)^4 + log(p1+4)^2 + (p2-pi)^2';
    case 3
        symb.Ob = '(p1-1)^4 + (p2-3)^2 + p1^2*p2';
end

symb.o.ny = 1; % number of states
symb.o.np = 2; % number of parameters
symb.o.output = 2; % interp1 compatible

p.t0 = 0; p.tf = 1;

% simple bounds
LB(1).right = 3; LB(1).matrix = [1e-10,1e-10]; % parameters

% combine structures
setup.symb = symb; setup.LB = LB;
setup.t0 = p.t0; setup.tf = p.tf; setup.p = p;

opts = [];
opts.general.displevel = 1;
opts.general.plotflag = 1;
opts.dt.nt = 2;

% solve
t1 = tic;
[T,U,Y,P1,F1,in,opts] = DTQP_solve(setup,opts);
toc(t1)

disp("F1"); disp(F1); disp("P1"); disp(P1(:)')

% fmincon options
options = optimoptions('fmincon','Algorithm','sqp','Display','iter',...
    'ScaleProblem',false);

% solve with fmincon and an sqp algorithm
switch testnum
    case 1
        [P2,F2] = fmincon(@(x) (x(1)-1)^4 + x(1)^2 + (x(2)-pi)^2,[1,1],[],[],[],[],[],[],[],options);
    case 2
        [P2,F2] = fmincon(@(x) (x(1)-1)^4 + log(x(1)+4)^2 + (x(2)-pi)^2,[1,1],[],[],[],[],[],[],[],options);
    case 3
        [P2,F2] = fmincon(@(x) (x(1)-1)^4 + (x(2)-3)^2 + x(1)^2*x(2),[1,1],[],[],[],[],[],[],[],options);
end

disp("F2"); disp(F2); disp("P2"); disp(P2(:)')