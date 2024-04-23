%--------------------------------------------------------------------------
% LQRScalarTransfer.m
%
% Similar behavior to the "Hyper-Sensitive" problem:
% A. V. Rao and K. D. Mease, Eigenvector Approximate Dichotomic Basis
% Methods for Solving Hyper-Sensitive Optimal Control Problems, Optimal
% Control Applications and Methods, Vol. 21, No. 1., January-February,
% 2000, pp. 1-17.
%
% "Energy-Optimal Control" problem in reference below is a special case:
% H. P. Geering, Optimal Control with Engineering Applications, Springer,
% 2007, pp. 46-48, doi: 10.1007/978-3-540-69438-0
%
% "Mass-Spring" problem included in GPOPS-II is a special case.
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

% auxiliary data:

% "Hyper-Sensitive"-like behavior
tf = 10000; % final time, requires high-precision solution
auxdata.a = -1; % state matrix
auxdata.b = 1; % input matrix
auxdata.c = 1.5; % initial state
auxdata.d = 1; % final state
auxdata.q = 1; % quadratic state cost
auxdata.r = 1; % quadratic control cost
setup.auxdata = auxdata;

% "Energy-Optimal Control"
% tf = 5; % final time, tf > 0
% auxdata.a = 2; % state matrix, a > 0
% auxdata.b = 3; % input matrix, b > 0
% auxdata.c = 1; % initial state, c > 0
% auxdata.d = 5; % final state, d < exp(a*tf)*c
% auxdata.q = 0; % quadratic state cost
% auxdata.r = 1/2; % quadratic control cost
% setup.auxdata = auxdata;

% "Mass-Spring" Problem
% tf = pi/2; % final time
% auxdata.a = 0; % state matrix
% auxdata.b = 1; % input matrix
% auxdata.c = 0; % initial state
% auxdata.d = 1; % final state
% auxdata.q = -1; % quadratic state cost
% auxdata.r = 1; % quadratic control cost
% setup.auxdata = auxdata;

% time horizon
setup.t0 = 0; setup.tf = tf;

% system dynamics
setup.lq.dynamics.A = auxdata.a;
setup.lq.dynamics.B = auxdata.b;

% Lagrange term
L(1).left = 1; % control variables
L(1).right = 1; % control variables
L(1).matrix = auxdata.r;
L(2).left = 2; % state variables
L(2).right = 2; % state variables
L(2).matrix = auxdata.q;
setup.lq.lagrange = L;

% simple bounds
LB(1).right = 4;  LB(1).matrix = auxdata.c; % initial states
UB(1).right = 4;  UB(1).matrix = auxdata.c; % initial states
LB(2).right = 5;  LB(2).matrix = auxdata.d;  % final states
UB(2).right = 5;  UB(2).matrix = auxdata.d; % final states
setup.lq.ub = UB; setup.lq.lb = LB;

end


%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,setup,opts)

% outputs
[O,sol] = LQRScalarTransfer_output(T,U,Y,P,F,in,opts);

% plots
LQRScalarTransfer_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts
% test number
num = 2;

switch num
case 1
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'CGL';
    opts.dt.nt = 2000; % number of nodes
case 2
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 10; % number of nodes
    opts.solver.tolerance = 1e-15;
    opts.solver.display = 'none';
    opts.dt.meshr.method = 'SS-BETTS';
    opts.dt.meshr.tolerance = 1e-6;
    opts.dt.meshr.ntmaxinterval = 5;
end

end