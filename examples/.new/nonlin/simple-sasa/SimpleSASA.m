%--------------------------------------------------------------------------
% SimpleSASA.m
% D. R. Herber and J. T. Allison, "Unified Scaling of Dynamic Optimization
% Design Formulations," in Volume 2A: 43rd Design Automation Conference,
% 2017, doi: 10.1115/detc2017-67676
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
auxdata.umax = 1;
auxdata.J = 1;
auxdata.t0 = 0;
auxdata.tf = 2;
setup.auxdata = auxdata;

% time horizon
setup.t0 = auxdata.t0; setup.tf = auxdata.tf;

% number of controls, states, and parameters
setup.counts.nu = 1;
setup.counts.nx = 2;
setup.counts.np = 1;

% Mayer term
M(1).right = 5; M(1).left = 0; M(1).matrix = [-1,0]; % final states
setup.lq.mayer = M;

% system dynamics
setup.nonlin.dynamics = '[y2; -p1/J*y1 + u1/J]';

% symbolic data for nonlin
setup.nonlin.data.symbols = 'J';
setup.nonlin.data.values = [auxdata.J];

% simple bounds
UB(1).right = 4; UB(1).matrix = [0,0]; % initial states
LB(1).right = 4; LB(1).matrix = [0,0];
UB(2).right = 5; UB(2).matrix = [inf 0]; % final states
LB(2).right = 5; LB(2).matrix = [-inf 0];
UB(3).right = 1; UB(3).matrix = auxdata.umax; % controls
LB(3).right = 1; LB(3).matrix = -auxdata.umax;
setup.lq.ub = UB; setup.lq.lb = LB;

% guess
Y0 = [[0,0];[1 0]];
U0 = [0;0];
P0 = [1;1];
setup.method.guess.X = [U0,Y0,P0];

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = SimpleSASA_output(T,U,Y,P,F,in,opts);

% plots
SimpleSASA_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts

% test number
num = 1;

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
end

end