%--------------------------------------------------------------------------
% BrysonHo59.m
% pp. 59-63 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
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
auxdata.t0 = 0; 
auxdata.tf = 2.05;
auxdata.a = 1;
auxdata.h = 1;
setup.auxdata = auxdata;

% time horizon
setup.t0 = auxdata.t0; setup.tf = auxdata.tf;

% number of controls, states, and parameters
counts.nu = 1;
counts.nx = 4;
counts.np = 0;
counts.nv = 0;
setup.counts = counts;

% Mayer term
M(1).left = 0; M(1).right = 5; M(1).matrix = [-1,0,0,0];
setup.lq.mayer = M;

% system dynamics
setup.nonlin.dynamics = '[a*cos(u1); a*sin(u1); y1; y2]';
setup.nonlin.data.symbols = 'a';
setup.nonlin.data.values = [auxdata.a];

% simple bounds
UB(1).right = 4; UB(1).matrix = [0,0,0,0]; % initial states
LB(1).right = 4; LB(1).matrix = [0,0,0,0];
UB(2).right = 5; UB(2).matrix = [inf,0,inf,auxdata.h]; % final states
LB(2).right = 5; LB(2).matrix = [-inf,0,-inf,auxdata.h];
UB(3).right = 1; UB(3).matrix = pi; % controls
LB(3).right = 1; LB(3).matrix = -pi;
UB(4).right = 2; UB(4).matrix = [inf,inf,inf,inf]; % states
LB(4).right = 2; LB(4).matrix = [-inf,-inf,-inf,0];
setup.lq.ub = UB; setup.lq.lb = LB;

% guess
Y0 = [[0,0,0,0];[0,0,0,auxdata.h]];
U0 = [0;0];
setup.method.guess.X = [U0,Y0];

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = BrysonHo59_output(T,U,Y,P,F,in,opts);

% plots
BrysonHo59_plot(T,U,Y,P,F,in,opts,sol)

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
    opts.solver.tolerance = 1e-8;
    opts.solver.maxiters = 30;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 60; % number of nodes
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
case 3 % qlin method
    opts.general.plotflag = 1; % create the plots
    opts.general.displevel = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100; % number of nodes
    opts.solver.maxiters = 200;
    opts.solver.display = 'none';
    opts.method.form = 'qlin';
    opts.method.trustregionflag = true;
    opts.method.sqpflag = false;
    opts.method.delta = 1e9;
    opts.method.improveguess = true;
end

end