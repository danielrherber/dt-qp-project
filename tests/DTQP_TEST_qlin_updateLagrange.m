%--------------------------------------------------------------------------
% DTQP_TEST_qlin_updateLagrange.m
% Test function for DTQP_QLIN_update_lagrange.m
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
clc; clear; close all

% tests = 1:5;
tests = 3;

% go through the tests
for k = 1:length(tests)

    % test setup
    switch tests(k)
        %------------------------------------------------------------------
        case 1
        symb.Ob = '1';
        n.nu = 2;
        n.ny = 3;
        n.np = 1;
        %------------------------------------------------------------------
        case 2
        symb.Ob = 'y1*u1';
        n.nu = 1;
        n.ny = 1;
        n.np = 0;
        %------------------------------------------------------------------
        case 3
        symb.Ob = 'y1^2 - y2^2 - u1^2 + u2^2 + y1*u1';
        n.nu = 2;
        n.ny = 2;
        n.np = 0;
        %------------------------------------------------------------------
        case 4
        symb.Ob = 'y1*u1^2 ';
        n.nu = 1;
        n.ny = 2;
        n.np = 0;
        %------------------------------------------------------------------
        case 5
        symb.Ob = 'y1^3 + y1*y2 + p1*u1';
        n.nu = 2;
        n.ny = 3;
        n.np = 1;
    end

    % problem structure
    [setup,opts,T,X,param] = problem(symb,n);
    symb = setup.symb; L = setup.symb.L; o = n;

    % run the test and time
    setup = DTQP_QLIN_update_lagrange(setup,L.H,L.G,L.C,o,T,X,param);

    % test analysis
    disp(strcat("-> Case: ",string(tests(k))))
    for i = 1:length(setup.L)
        disp(strcat("Left: ",string(setup.L(i).left)))
        disp(strcat("Right: ",string(setup.L(i).right)))
        disp(setup.L(i).matrix)
        t = linspace(T(1),T(end),3)';
        A = DTQP_tmultiprod(setup.L(i).matrix,[],t);
        disp(' ')
    end

end

% problem structure
function [setup,opts,T,X,param] = problem(symb,n)

n.output = 2;
setup.symb = symb;
setup.n = n;
setup.t0 = 0;
setup.tf = 1;
param = [];

% options
opts = [];
opts.general.displevel = false;
opts.dt.nt = 4;
opts.method.form = 'qlin';

% initialize some stuff
[setup,opts] = DTQP_default_opts(setup,opts);

% run the test and time
setup = DTQP_QLIN_initialize(setup,opts); % FIX

% initial guess
[T,U,Y,P] = DTQP_QLIN_guess(setup,opts,n);

% construct previous solution vector
P = repelem(P',opts.dt.nt,1);
X = [U,Y,P];

end