%--------------------------------------------------------------------------
% BrysonHo166.m
% pp. 166-167 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
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
auxdata.x0 = -0.5; auxdata.v0 = 1; % other
setup.auxdata = auxdata;

% time horizon
setup.t0 = 0; setup.tf = 20;

% system dynamics
setup.lq.dynamics.A = [0 1;-1 0];
setup.lq.dynamics.B = [0;1];

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = 1/2; % 1/2*u^2

% simple bounds
LB(1).right = 4;  LB(1).matrix = [auxdata.x0;auxdata.v0]; % initial states
UB(1).right = 4;  UB(1).matrix = [auxdata.x0;auxdata.v0];
LB(2).right = 5;  LB(2).matrix = [0;0]; % final states
UB(2).right = 5;  UB(2).matrix = [0;0];
setup.lq.ub = UB; setup.lq.lb = LB;

end


%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = BrysonHo166_output(T,U,Y,P,F,in,opts);

% plots
BrysonHo166_plot(T,U,Y,P,F,in,opts,sol)

end



%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts
% test number
num = 2;

switch num
case 1
    % default parameters
    opts.general.plotflag = 1; % create the plots
    opts.general.saveflag = 0;
    opts.general.displevel = 2;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100; % number of nodes
    opts.method.reordervariables = 0;
    opts.solver.function = 'built-in';
    opts.solver.tolerance = 1e-12;
    opts.solver.maxiters = 200;
case 2
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 5;
    opts.method.reordervariables = 0;
    opts.solver.function = 'built-in';
    opts.solver.tolerance = 1e-15;
    opts.solver.maxiters = 200;
    opts.solver.display = 'none';
    opts.dt.meshr.method = 'SS-BETTS';
    opts.dt.meshr.tolerance = 1e-4;
case 3
    opts.dt.defects = 'PS-MI';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 10; % number of intervals
    opts.dt.nn = 12; % polynomial order in each interval
end

end