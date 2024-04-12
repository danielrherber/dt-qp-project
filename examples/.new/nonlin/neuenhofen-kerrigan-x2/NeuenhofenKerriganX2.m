%--------------------------------------------------------------------------
% NeuenhofenKerriganX2.m
% X2 from M. P. Neuenhofen and E. C. Kerrigan, "An integral penalty-barrier
% direct transcription method for optimal control," 2020, arXiv:2009.06217
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
setup.t0 = 0; setup.tf = 1;

% number of controls, states, and parameters
setup.counts.nu = 1;
setup.counts.nx = 2;

% Mayer term
M(1).left = 0; M(1).right = 5; M(1).matrix = [0 1]; % states
setup.lq.mayer = M;

% system dynamics
setup.nonlin.dynamics = '[u1/(2*y1); 4*y1^4 + u1^2]';

% simple bounds
UB(1).right = 4; UB(1).matrix = [1 0]; % initial states
LB(1).right = 4; LB(1).matrix = [1 0];
LB(2).right = 1; LB(2).matrix = -1; % controls
LB(3).right = 2; LB(3).matrix = [sqrt(0.4),-inf]; % states
setup.lq.ub = UB; setup.lq.lb = LB;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = NeuenhofenKerriganX2_output(T,U,Y,P,F,in,opts);

% plots
NeuenhofenKerriganX2_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts

% test number
num = 1;

switch num
case 1
    opts.general.plotflag = 1; % create the plots
    opts.general.saveflag = false;
    opts.general.displevel = 2;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 200; % number of nodes
    opts.solver.tolerance = 1e-12;
    opts.solver.maxiters = 200;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
case 2
    opts.general.plotflag = 1; % create the plots
    opts.general.saveflag = false;
    opts.general.displevel = 2;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 50; % number of nodes
    opts.solver.tolerance = 1e-12;
    opts.solver.maxiters = 200;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
end

end