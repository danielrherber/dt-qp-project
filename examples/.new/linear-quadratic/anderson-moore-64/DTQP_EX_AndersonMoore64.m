%--------------------------------------------------------------------------
% AndersonMoore64.m
% p. 64 of B. D. O. Anderson and J. B. Moore, Optimal Control: Linear
% Quadratic Methods. Prentice-Hall, 1989, isbn: 0136386512
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

% tunable parameters
auxdata.x0 = 5; % initial state
setup.auxdata = auxdata;

% time horizon
setup.t0 = 0; setup.tf = 10;

% Lagrange term
L(1).left = 1; % controls
L(1).right = 1; % controls
L(1).matrix = {@(t) 2*exp(-t)};
L(2).left = 2; % states
L(2).right = 2; % states
L(2).matrix = {@(t) 0.5*exp(-t)};
setup.lq.lagrange = L;

% system dynamics
setup.lq.dynamics.A = 1/2; setup.lq.dynamics.B = 1;

% simple bounds
LB(1).right = 4; % initial states
LB(1).matrix = auxdata.x0;
UB(1).right = 4; % initial states
UB(1).matrix = auxdata.x0;
setup.lq.ub = UB; setup.lq.lb = LB;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = AndersonMoore64_output(T,U,Y,P,F,in,opts);

% plots
AndersonMoore64_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts

% test number
num = 1;

switch num
    case 1
        opts.dt.defects = 'PS';
        opts.dt.quadrature = 'G';
        opts.dt.mesh = 'LGL';
        opts.dt.nt = 15;
    case 2
        opts.dt.quadrature = 'CQHS';
        opts.dt.defects = 'HS';
        opts.dt.mesh = 'ED';
        opts.dt.nt = 100;
    case 3
        opts.dt.quadrature = 'CTR';
        opts.dt.defects = 'TR';
        opts.dt.mesh = 'ED';
        opts.dt.nt = 100;
end

end