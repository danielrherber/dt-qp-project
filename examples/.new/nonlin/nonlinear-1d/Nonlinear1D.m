%--------------------------------------------------------------------------
% Nonlinear1D.m
% D. Garg et al., "Direct Trajectory Optimization and Costate Estimation of
% General Optimal Control Problems Using a Radau Pseudospectral Method," in
% AIAA Guidance, Navigation, and Control Conference, 2009,
% doi: 10.2514/6.2009-5989
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
setup.t0 = 0; setup.tf = 5;

% number of controls, states, and parameters
setup.counts.nu = 1;
setup.counts.nx = 1;

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = 1/2; % u^2
L(2).left = 0; L(2).right = 2; L(2).matrix = 1/2; % x
setup.lq.lagrange = L;

% system dynamics
setup.nonlin.dynamics = '[2*y1 + 2*u1*sqrt(y1)]';

% simple bounds
UB(1).right = 4; UB(1).matrix = 2; % initial states
LB(1).right = 4; LB(1).matrix = 2;
UB(2).right = 5; UB(2).matrix = 1; % final states
LB(2).right = 5; LB(2).matrix = 1;
LB(3).right = 2; LB(3).matrix = 0; % optional constraint for sqrt(y1)
setup.lq.ub = UB; setup.lq.lb = LB;

% guess
Y0 = [2;1];
U0 = [0;0];
setup.method.guess.X = [U0,Y0];

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = Nonlinear1D_output(T,U,Y,P,F,in,opts);

% plots
Nonlinear1D_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts

% test number
num = 2;

switch num
case 1
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100; % number of nodes
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.solver.tolerance = 1e-10;
    opts.dt.meshr.method = 'RICHARDSON-DOUBLING';
    opts.dt.meshr.tolerance = 1e-6;
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 8; % number of nodes
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.solver.tolerance = 1e-10;
    opts.dt.meshr.method = 'RICHARDSON-DOUBLING';
    opts.dt.meshr.tolerance = 1e-6;
end

end