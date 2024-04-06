%--------------------------------------------------------------------------
% Rayleigh.m
% H. Maurer and D. Augustin, "Sensitivity Analysis and Real-Time Control of
% Parametric Optimal Control Problems Using Boundary Value Methods," Online
% Optimization of Large Scale Systems. Springer Berlin Heidelberg, pp.
% 17–55, 2001. doi: 10.1007/978-3-662-04331-8_2
%
% Also, Example 4.6 in J. T. Betts, "Practical Methods for Optimal Control
% and Estimation Using Nonlinear Programming." Society for Industrial and
% Applied Mathematics, Jan. 01, 2010. doi: 10.1137/1.9780898718577
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
auxdata.p = 0.14;
auxdata.t0 = 0;
auxdata.tf = 4.5;

% time horizon
setup.t0 = auxdata.t0; setup.tf = auxdata.tf;

% number of controls, states, and parameters
counts.nu = 1;
counts.nx = 2;
counts.np = 0;
counts.nv = 0;
setup.counts = counts;

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = 1; % u^2
L(2).left = 2; L(2).right = 2; L(2).matrix = [1,0;0,0]; % y1^2
setup.lq.lagrange = L;

% system dynamics
setup.nonlin.dynamics = '[y2; -y1+y2*(1.4-p*y2^2)+4*u1]';
setup.nonlin.data.symbols = 'p';
setup.nonlin.data.values = [auxdata.p];

% simple bounds
UB(1).right = 4; UB(1).matrix = [-5,-5]; % initial states
LB(1).right = 4; LB(1).matrix = [-5,-5];
UB(2).right = 5; UB(2).matrix = [0,0]; % final states
LB(2).right = 5; LB(2).matrix = [0,0];
UB(3).right = 1; UB(3).matrix = 1; % controls
LB(3).right = 1; LB(3).matrix = -1;
setup.lq.ub = UB; setup.lq.lb = LB;

% guess
Y0 = [[-5,-5];[0,0]];
U0 = [0;0];
setup.method.guess.X = [U0,Y0];

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = Rayleigh_output(T,U,Y,P,F,in,opts);

% plots
Rayleigh_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts

% test number
num = 1;

switch num
case 1
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100;
    opts.solver.tolerance = 1e-12;
    opts.solver.display = 'iter';
    opts.method.form = 'nonlinearprogram';
case 2
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 50;
    opts.solver.tolerance = 1e-12;
    opts.solver.display = 'iter';
    opts.method.form = 'nonlinearprogram';
end

end