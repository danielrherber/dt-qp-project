%--------------------------------------------------------------------------
% TavallaeiTousi1.m
% Based on the example in:
% M. A. Tavallaei and B. Tousi, "Closed Form Solution to an Optimal Control
% Problem by Orthogonal Polynomial Expansion," American J. of Engineering
% and Applied Sciences, vol. 1, no. 2, pp. 104-109, 2008.
% doi: 10.3844/ajeassp.2008.104.109
%--------------------------------------------------------------------------
% NOTE: the results presented in the paper don't seem to be for the
% matrices in the numerical example. Instead
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
auxdata.y0 = [-4;4];
auxdata.k = 3;
auxdata.r = 0.5;
auxdata.Q = [1,0;0,2];
auxdata.A = [2,0;1,-1];
auxdata.B = [1;0]; % modified from reference

% ns = 2; nu = 1;
% auxdata.tf = 5;
% auxdata.y0 = 5*(rand(ns,1)-0.5);
% auxdata.k = 2*pi;
% auxdata.r = 5*diag(rand(nu,1));
% auxdata.Q = 5*diag(rand(ns,1));
% auxdata.A = 5*(rand(ns)-0.5);
% auxdata.B = 5*(rand(ns,nu)-0.5);

setup.auxdata = auxdata;

% time horizon
setup.t0 = 0; setup.tf = 5;

% system dynamics
setup.lq.dynamics.A = auxdata.A;
setup.lq.dynamics.B = auxdata.B;

% Lagrange term
L(1).left = 1;  % control variables
L(1).right = 1; % control variables
L(1).matrix = auxdata.r;
L(2).left = 2;  % state variables
L(2).right = 2; % state variables
L(2).matrix = {'prod',{@(t) t.^auxdata.k},auxdata.Q};
setup.lq.lagrange = L;

% simple bounds
LB(1).right = 4;  LB(1).matrix = auxdata.y0; % initial states
UB(1).right = 4;  UB(1).matrix = auxdata.y0; % initial states
setup.lq.ub = UB; setup.lq.lb = LB;

end


%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,setup,opts)

% outputs
[O,sol] = TavallaeiTousi1_output(T,U,Y,P,F,in,opts,setup);

% plots
TavallaeiTousi1_plot(T,U,Y,P,F,in,opts,sol)

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
    opts.dt.nt = 20;
case 2
    opts.dt.quadrature = 'CQHS';
    opts.dt.defects = 'HS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100;
case 3
    opts.dt.quadrature = 'CTR';
    opts.dt.defects = 'TR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100;
end

end