%--------------------------------------------------------------------------
% BrysonHo64.m
% p. 64 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
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
auxdata.a = 1;
auxdata.l = 0.45;
auxdata.t0 = -auxdata.l;
auxdata.tf = auxdata.l;
setup.auxdata = auxdata;

% number of controls, states, and parameters
setup.counts.nu = 1;
setup.counts.nx = 1;

% time horizon
setup.t0 = auxdata.t0; setup.tf = auxdata.tf;

% objective function
setup.nonlin.lagrange = '2*pi*y1*sqrt(1+u1^2)';

% dynamics
setup.lq.dynamics.A = 0;
setup.lq.dynamics.B = 1;

% bounds
UB(1).right = 4; UB(1).matrix = auxdata.a;
LB(1).right = 4; LB(1).matrix = auxdata.a;
UB(2).right = 5; UB(2).matrix = auxdata.a;
LB(2).right = 5; LB(2).matrix = auxdata.a;
LB(3).right = 2; LB(3).matrix = 0;
setup.lq.ub = UB; setup.lq.lb = LB;

% guess
Y0 = [0;0];
U0 = [0;0];
setup.method.guess.X = [U0,Y0];

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = BrysonHo64_output(T,U,Y,P,F,in,opts);

% plots
BrysonHo64_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts

% test number
num = 2;

switch num
case 1
    opts.general.plotflag = 1; % create the plots
    opts.general.saveflag = false;
    opts.general.displevel = 2;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 30; % number of nodes
    opts.solver.tolerance = 1e-8;
    opts.solver.maxiters = 200;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 10; % number of nodes
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
end

end