%--------------------------------------------------------------------------
% SimpleCoDesignTransfer.m
% D. R. Herber, J. T. Allison. 'Nested and simultaneous solution strategies
% for general combined plant and control design problems.' ASME Journal of
% Mechanical Design, 141(1), p. 011402, Jan 2019. doi: 10.1115/1.4040705
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
setup = createSetup(opts);

% solve with DTQP
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

% post-processing
[O,sol] = postProcessing(T,U,Y,P,F,in,opts);

%--------------------------------------------------------------------------
% create setup structure
%--------------------------------------------------------------------------
function setup = createSetup(opts)

% initialize struct
setup = DTQP_setup_initialize;

% auxiliary data
auxdata.x0 = 1; auxdata.v0 = 2;
setup.auxdata = auxdata;

% time horizon
setup.t0 = 0;  setup.tf = 1;

% number of controls, states, and parameters
setup.counts.nu = 1;
setup.counts.nx = 2;
setup.counts.np = 1;

% Lagrange term
if ~isfield(opts.method,'olqflag') || opts.method.olqflag
    L(1).left = 1; L(1).right = 1; L(1).matrix = 1; % u^2
    setup.lq.lagrange = L;
else
    setup.nonlin.lagrange = 'u1^2';
end

% system dynamics
setup.nonlin.dynamics = '[y2;-p1*y1 + u1]';

% simple bounds
UB(1).right = 4; UB(1).matrix = [auxdata.x0;auxdata.v0]; % initial states
LB(1).right = 4; LB(1).matrix = [auxdata.x0;auxdata.v0];
UB(2).right = 5; UB(2).matrix = [0;0]; % final states
LB(2).right = 5; LB(2).matrix = [0;0];
setup.lq.ub = UB; setup.lq.lb = LB;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = SimpleCoDesignTransfer_output(T,U,Y,P,F,in,opts);

% plots
SimpleCoDesignTransfer_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts

% test number
num = 3;

switch num
case 1
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 200; % number of nodes
    opts.method.form = 'nonlinearprogram';
    opts.solver.function = 'ipfmincon';
    opts.solver.display = 'iter';
    opts.method.olqflag = true;
    opts.method.derivativeflag = true;
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 20; % number of nodes
    opts.method.form = 'nonlinearprogram';
    opts.solver.function = 'ipfmincon';
    opts.solver.display = 'iter';
    opts.method.olqflag = true;
    opts.method.derivativeflag = true;
case 3
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.method.maxiters = 5000;
    opts.dt.nt = 200; % number of nodes
    opts.method.form = 'qlin';
    opts.solver.display = 'none';
    opts.method.trustregionflag = false;
    opts.method.improveguess = false;
end

end