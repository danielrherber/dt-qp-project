%--------------------------------------------------------------------------
% OptimalConsumption.m
% L. C. Evans, *An Introduction to Mathematical Optimal Control Theory*.
% Department of Mathematics, University of California, Berkeley,
% Version 0.2, pp. 51-52
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
setup.t0 = 0; setup.tf = 1.5;

% number of controls, states, and parameters
setup.counts.nu = 1;
setup.counts.nx = 1;

% Lagrange term
L(1).left = 0; L(1).right = 2; L(1).matrix = -1;
L(2).left = 1; L(2).right = 2; L(2).matrix = 1;
setup.lq.lagrange = L;

% system dynamics
setup.nonlin.dynamics = '[y1*u1]';

% simple bounds
UB(1).right = 1; UB(1).matrix = 1; % controls
LB(1).right = 1; LB(1).matrix = 0;
LB(2).right = 4; LB(2).matrix = 1; % initial states
UB(2).right = 4; UB(2).matrix = 1;
setup.lq.ub = UB; setup.lq.lb = LB;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = OptimalConsumption_output(T,U,Y,P,F,in,opts);

% plots
OptimalConsumption_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts

% test number
num = 1;

switch num
case 1 % ipfmincon method
    opts.general.plotflag = 1; % create the plots
    opts.general.saveflag = false;
    opts.general.displevel = 2;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100; % number of nodes
    opts.solver.tolerance = 1e-15;
    opts.solver.maxiters = 200;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
case 2 % qlin method
    opts.general.plotflag = 1; % create the plots
    opts.general.displevel = 2;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100; % number of nodes
    opts.solver.maxiters = 200;
    opts.solver.display = 'none';
    opts.method.form = 'qlin';
    opts.method.trustregion = false;
    opts.method.sqpflag = false;
end


end