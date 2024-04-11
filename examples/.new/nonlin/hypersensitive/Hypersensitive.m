%--------------------------------------------------------------------------
% Hypersensitive.m
% A. V. Rao and K. D. Mease, "Eigenvector approximate dichotomic basis
% method for solving hyper-sensitive optimal control problems", *Optimal
% Control Applications and Methods*, vol. 21, no. 1, pp. 1-19, 2000,
% doi: 10.1002/(SICI)1099-1514(200001/02)21:1<1::AID-OCA646>3.0.CO;2-V
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

% time horizon
setup.t0 = 0; setup.tf = 10000;

% number of controls, states, and parameters
setup.counts.nu = 1;
setup.counts.nx = 1;

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = 1/2;
L(2).left = 2; L(2).right = 2; L(2).matrix = 1/2;
setup.lq.lagrange = L;

% system dynamics
setup.nonlin.dynamics = '[-y1^3+u1]';

% initial and final state value constraints
UB(1).right = 4; UB(1).matrix = 1.5; % initial states
LB(1).right = 4; LB(1).matrix = 1.5;
UB(2).right = 5; UB(2).matrix = 1; % final states
LB(2).right = 5; LB(2).matrix = 1;
setup.lq.ub = UB; setup.lq.lb = LB;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = Hypersensitive_output(T,U,Y,P,F,in,opts);

% plots
Hypersensitive_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts

% test number
num = 1;

switch num
case 1 % ipfmincon method
    opts.general.displevel = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    S2 = 10; N2 = 300; p.t0 = 0; p.tf = 10000;
    opts.dt.t = unique([linspace(p.t0,p.t0+S2,N2),linspace(p.t0+S2,p.tf-S2,N2),linspace(p.tf-S2,p.tf,N2)]);
    opts.dt.mesh = 'USER';
    opts.dt.nt = length(opts.dt.t);
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
case 2 % qlin method
    opts.general.displevel = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    S2 = 10; N2 = 300; p.t0 = 0; p.tf = 10000;
    opts.dt.t = unique([linspace(p.t0,p.t0+S2,N2),linspace(p.t0+S2,p.tf-S2,N2),linspace(p.tf-S2,p.tf,N2)]);
    opts.dt.mesh = 'USER';
    opts.dt.nt = length(opts.dt.t);
    opts.solver.display = 'none';
    opts.method.form = 'qlin';
    opts.method.trustregion = false;
    opts.method.sqpflag = false;
end

end