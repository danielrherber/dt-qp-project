%--------------------------------------------------------------------------
% BrysonHo248.m
% pp. 248-250 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
%--------------------------------------------------------------------------
% This is a singular optimal control problem
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
auxdata.alpha = 1; auxdata.beta = 1; auxdata.gamma = 20;
setup.auxdata = auxdata;

% time horizon
setup.t0 = 0; setup.tf = 1;

% system dynamics
setup.lq.dynamics.A = [0 1; 0 0];
setup.lq.dynamics.B = [1;-1];

% Lagrange term
L(1).left = 2; L(1).right = 2; L(1).matrix = [1/2,0;0,0]; % state variables
setup.lq.lagrange = L;

% simple bounds
LB(1).right = 4; LB(1).matrix = [auxdata.alpha;auxdata.beta]; % initial states
UB(1).right = 4; UB(1).matrix = [auxdata.alpha;auxdata.beta]; % initial states
LB(2).right = 5; LB(2).matrix = [0;0]; % final states
UB(2).right = 5; UB(2).matrix = [0;0]; % final states
LB(3).right = 1; LB(3).matrix = -auxdata.gamma; % control
UB(3).right = 1; UB(3).matrix = auxdata.gamma; % control
setup.lq.ub = UB; setup.lq.lb = LB;

end


%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,setup,opts)

% outputs
[O,sol] = BrysonHo248_output(T,U,Y,P,F,in,opts);

% plots
BrysonHo248_plot(T,U,Y,P,F,in,opts,sol)

end




%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts

% test number
num = 1;

switch num
case 1
    opts.dt.quadrature = 'CEF';
    opts.dt.defects = 'ZO';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100;
case 2
    opts.dt.quadrature = 'G';
    opts.dt.defects = 'PS';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 100;
case 3
    opts.dt.quadrature = 'CQHS';
    opts.dt.defects = 'HS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100;
case 4
    opts.dt.quadrature = 'CQHS';
    opts.dt.defects = 'HS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 10;
    opts.solver.display = 'none';
    opts.solver.tolerance = 1e-15;
    opts.dt.meshr.method = 'SS-BETTS';
    opts.dt.meshr.tolerance = 1e-2;
end

end