%--------------------------------------------------------------------------
% BrysonHo116.m
% pp. 116-117 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
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
auxdata.x0 = 0.1; auxdata.v0 = 0.5;
setup.auxdata = auxdata;

% time horizon
setup.t0 = 0; setup.tf = 1.5;

% system dynamics
setup.lq.dynamics.A = [0 1; 0 0];
setup.lq.dynamics.B = [0 0; 1 0];

% Lagrange term
L(1).right = 1;  L(1).matrix = [0,1];
setup.lq.lagrange = L;

% simple bounds
LB(1).right = 4; LB(1).matrix = [auxdata.x0;auxdata.v0]; % initial state
UB(1).right = 4; UB(1).matrix = [auxdata.x0;auxdata.v0];
LB(2).right = 5; LB(2).matrix = [0;0]; % final state
UB(2).right = 5; UB(2).matrix = [0;0];
LB(3).right = 1; LB(3).matrix = [-1;-Inf];% control
UB(3).right = 1; UB(3).matrix = [1;Inf]; 
setup.lq.ub = UB; setup.lq.lb = LB;


% absolute value approximation
inequality(1).linear(1).right = 1;
inequality(1).linear(1).matrix = [-1; -1];
inequality(1).b = 0;
inequality(2).linear(1).right = 1;
inequality(2).linear(1).matrix = [1;-1];
inequality(2).b = 0;
setup.lq.inequality = inequality;

end 

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = BrysonHo116_output(T,U,Y,P,F,in,opts);

% plots
BrysonHo116_plot(T,U,Y,P,F,in,opts,sol)

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
    opts = [];
case 2
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 5;
    opts.solver.tolerance = 1e-15;
    opts.solver.display = 'none';
    opts.dt.meshr.method = 'SS-BETTS';
    opts.dt.meshr.tolerance = 1e-7;
end

end