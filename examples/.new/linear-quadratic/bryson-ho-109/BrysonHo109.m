%--------------------------------------------------------------------------
% BrysonHo109.m
% pp. 109-110 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
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
[O,sol] = postProcessing(T,U,Y,P,F,in,setup,opts);

%--------------------------------------------------------------------------
% create setup structure
%--------------------------------------------------------------------------
function setup = createSetup

% initialize struct
setup = DTQP_setup_initialize;

% auxiliary data
auxdata.x0 = 1; auxdata.a = 2; % 1
setup.auxdata = auxdata;

% time horizon
setup.t0 = 0; setup.tf = 1;

% system dynamics
setup.lq.dynamics.A = 0;
setup.lq.dynamics.B{1,1} = @BrysonHo109_g;

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = 1/2; % 1/2*u^2
setup.lq.lagrange = L;

% Mayer term
M(1).left = 5; M(1).right = 5; M(1).matrix = auxdata.a^2/2; % a^2/2*xf^2
setup.lq.mayer = M;

% simple bounds
UB(1).right = 4; UB(1).matrix = auxdata.x0; % initial state
LB(1).right = 4; LB(1).matrix = auxdata.x0;
UB(2).right = 1; UB(2).matrix = 1; % control
LB(2).right = 1; LB(2).matrix = -1;
setup.lq.ub = UB; setup.lq.lb = LB;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,setup,opts)

% outputs
[O,sol] = BrysonHo109_output(T,U,Y,P,F,in,setup,opts);

% plots
BrysonHo109_plot(T,U,Y,P,F,in,opts,sol)

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
    opts.general.saveflag = 0;
    opts.general.displevel = 2;
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100; % number of nodes
    opts.method.reordervariables = 0;
    opts.solver.function = 'built-in';
    opts.solver.tolerance = 1e-15;
    opts.solver.maxiters = 200;
    opts.solver.display = 'iter';
case 2
    opts.general.plotflag = 1; % create the plots
    opts.general.saveflag = 0;
    opts.general.displevel = 2;
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 10; % number of nodes
    opts.method.reordervariables = 0;
    opts.solver.function = 'built-in';
    opts.solver.tolerance = 1e-15;
    opts.solver.maxiters = 200;
    opts.solver.display = 'none';
    opts.dt.meshr.method = 'SS-BETTS';
    opts.dt.meshr.tolerance = 1e-7;
end

end