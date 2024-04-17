%--------------------------------------------------------------------------
% DTQP2.m
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

% auxiliary data
auxdata.x0 = 1; auxdata.r = 1; auxdata.m = 1; 
auxdata.a = 1; auxdata.b = 1; auxdata.omega = pi;
setup.auxdata = auxdata;

% time horizon
setup.t0 = 0; setup.tf = 15;

% system dynamics
setup.lq.dynamics.A = auxdata.a;
setup.lq.dynamics.B = @(t) auxdata.b*sin(auxdata.omega*t);

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = auxdata.r; % control variables
setup.lq.lagrange = L;

% Mayer term
M(1).left = 5; M(1).right = 5; M(1).matrix = auxdata.m; % final states
setup.lq.mayer = M;

% simple bounds
LB(1).right = 4; LB(1).matrix = auxdata.x0; % initial states
UB(1).right = 4;  UB(1).matrix = auxdata.x0; % initial states
setup.lq.ub = UB; setup.lq.lb = LB;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = DTQP2_output(T,U,Y,P,F,in,opts);

% plots
DTQP2_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts
% test number
num = 2;

switch num
case 1
    opts.dt.nt = 1000; % number of time points
case 2
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 10;
    opts.solver.display = 'none';
    opts.dt.meshr.method = 'SS-BETTS';
    opts.dt.meshr.tolerance = 1e-5;
end

end