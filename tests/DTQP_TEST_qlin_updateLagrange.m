%--------------------------------------------------------------------------
% DTQP_TEST_qlin_updateLagrange.m
% Test function for DTQP_qlin_updateLagrange.m
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
clc; clear; close all

testnum = 3;

switch testnum
    case 1
        setup.symb.Ob = '1';

    case 2
        setup.symb.Ob = 'y1*u1';

    case 3
        setup.symb.Ob = 'y1^3 + y1*y2 + p1*u1';
end

o.nu = 2;
o.ny = 3;
o.np = 1;
setup.symb.o = o;
setup.p.t0 = 0;
setup.p.tf = 1;

param = [];

% options
opts = [];
opts.dt.nt = 10;

% initial guess
[T,U,Y,P] = DTQP_qlin_guess(setup,opts,o);

% construct previous solution vector
P = repelem(P',opts.dt.nt,1);
X = [U,Y,P];

% initialize some stuff
[setup,opts] = DTQP_default_opts(setup,opts);

% modify setup
setup = DTQP_qlin_initialize(setup,opts); % FIX
L = setup.symb.L;

setup = DTQP_qlin_updateLagrange(setup,L.H,L.G,L.C,o,T,X,param);

% display the L structure
for k = 1:length(setup.L)
    disp(setup.L(k))
end