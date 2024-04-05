%--------------------------------------------------------------------------
% BettsBiehnCampbell1.m
% J. T. Betts, N. Biehn, and S. L. Campbell, Convergence of Nonconvergent
% IRK Discretizations of Optimal Control Problems with State Inequality
% Constraints," SIAM Journal on Scientific Computing, vol. 23, no. 6,
% pp. 1981-2007, 2002. doi: 10.1137/S1064827500383044
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clearvars -except externalInput

% set options and standardize
if ~exist('externalInput','var')
    opts = localOpts;
end
DTQP_standardizedinputs2

% create setup structure
setup = createSetup;

% solve with DTQP
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

% post-processing
[O,sol] = postProcessing(T,U,Y,P,F,in,opts);

%--------------------------------------------------------------------------
% create setup structure
%--------------------------------------------------------------------------
function setup = createSetup

% initialize struct
setup = DTQP_setup_initialize;

% time horizon
t0 = 34/15; tf = 4;

% system dynamics
A = [0 1; 0 0]; B = [0;1];

% Lagrange term
L(1).left = 1; % controls
L(1).right = 1; % controls
L(1).matrix = 1e-3; % 1e-3*u^2
L(2).left = 2; % states
L(2).right = 2; % states
L(2).matrix = [1 0; 0 0]; % y1^2

% initial states
LB(1).right = 4; % initial states
LB(1).matrix = [302399/50625; 70304/3375];
UB(1).right = 4; % initial states
UB(1).matrix = [302399/50625; 70304/3375];

% simple bound path constraint
LB(2).right = 2; % states
LB(2).matrix = {@(t) 15 - (t-4).^4;-Inf};

% add
setup.auxdata = [];
setup.t0 = t0;
setup.tf = tf;
setup.lq.lagrange = L;
setup.lq.dynamics.A = A;
setup.lq.dynamics.B = B;
setup.lq.ub = UB;
setup.lq.lb = LB;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = BettsBiehnCampbell1_output(T,U,Y,P,F,in,opts);

% plots
BettsBiehnCampbell1_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts

% test number
num = 1;

switch num
    case 1
        opts.dt.defects = 'TR';
        opts.dt.quadrature = 'CTR';
        opts.dt.mesh = 'ED';
        opts.dt.nt = 110; % number of nodes
    case 2
        opts.dt.defects = 'PS';
        opts.dt.quadrature = 'G';
        opts.dt.mesh = 'LGL';
        opts.dt.nt = 11; % number of nodes
end

end