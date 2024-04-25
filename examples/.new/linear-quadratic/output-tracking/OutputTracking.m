%--------------------------------------------------------------------------
% OutputTracking.m
% p. 175 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
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

% tunable parameters
ny = 8; % number of states
nu = 3; % number of controls
no = 2; % number of outputs

% time horizon
setup.t0 = 0; setup.tf = 50;

rng(900911075,'twister') % specific random seed
x0 = ones(ny,1); % initial states
A = sprand(ny,ny,0.7,1); % state matrix
B = sprand(ny,nu,1,1); % input matrix
C = sprand(no,ny,1,1); % output matrix
R = 1e-2*eye(nu); % control penalty
Q = eye(no); % state penalty

% output to track
[o,W] = OutputTracking_o(no);

% auxiliary data
auxdata.x0 = x0; auxdata.A = A; auxdata.B = B; auxdata.C = C;
auxdata.R = R; auxdata.Q = Q; auxdata.W = W; auxdata.o = o;
setup.auxdata = auxdata;

% system dynamics
setup.lq.dynamics.A = A;
setup.lq.dynamics.B = B;

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = R; % u'*R*u
L(2).left = 2; L(2).right = 2; L(2).matrix = C'*Q*C; % (C*x)'*Q*(C*x)
L(3).left = 0; L(3).right = 2; L(3).matrix = {'prod',o',-2*Q,C}; % -2*o'*Q*C
L(4).left = 0; L(4).right = 0; L(4).matrix = {'prod',o',Q,o}; % o'*Q*o
setup.lq.lagrange = L;

% simple bounds
UB(1).right = 4; UB(1).matrix = x0; % initial states
LB(1).right = 4; LB(1).matrix = x0; % initial states
setup.lq.ub = UB; setup.lq.lb = LB;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,setup,opts)

% outputs
[O,sol] = OutputTracking_output(T,U,Y,P,F,in,setup,opts);

% plots
OutputTracking_plot(T,U,Y,P,F,in,opts,sol)

end



%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts
% test number
num = 1;

switch num
case 1
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 100; % number of nodes
case 2
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 1000; % number of nodes
case 3
    opts = [];
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 50;
    opts.solver.tolerance = 1e-15;
    opts.solver.display = 'none';
    opts.dt.meshr.method = 'SS-BETTS';
    opts.dt.meshr.tolerance = 1e-3;
end

end