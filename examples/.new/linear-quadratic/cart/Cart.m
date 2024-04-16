%--------------------------------------------------------------------------
% Cart.m
% pp. 124-126 of D. H. Ballard, An Introduction to Natural Computation. 
% MIT Press, 1999, isbn: 9780262522588
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

% time horizon
setup.t0 = 0; setup.tf = 1;

% system dynamics
setup.lq.dynamics.A = [0 1; 0 -1]; 
setup.lq.dynamics.B = [0;1];

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix(1,1) = 1/2; % 1/2*u.^2 control variables
setup.lq.lagrange = L;

% Mayer term
M(1).left = 0; M(1).right = 5; M(1).matrix = [-1,0];  % singleto-final states
setup.lq.mayer = M;

% simple bounds
LB(1).right = 4; LB(1).matrix = [0;0]; % initial states
UB(1).right = 4; UB(1).matrix = [0;0]; % initial states
setup.lq.ub = UB; setup.lq.lb = LB;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,setup,opts)

% outputs
[O,sol] = Cart_output(T,U,Y,P,F,in,opts);

% plots
Cart_plot(T,U,Y,P,F,in,opts,sol)

end





%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts
% test number
num = 1;

switch num
case 1
    opts.dt.nt = 100; % number of nodes
end

end