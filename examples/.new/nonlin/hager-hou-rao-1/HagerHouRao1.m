%--------------------------------------------------------------------------
% HagerHouRao1.m
% W. W. Hager, H. Hou, and A. V. Rao, "Convergence Rate for a Gauss
% Collocation Method Applied to Unconstrained Optimal Control," J Optim
% Theory Appl, vol. 169, no. 3, pp. 801–824, Mar. 2016,
% doi: 10.1007/s10957-016-0929-7
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
auxdata.t0 = 0; auxdata.tf = 1;
setup.auxdata = auxdata;

% time horizon
setup.t0 = 0; setup.tf = 1;

% number of controls, states, and parameters
setup.counts.nu = 2;
setup.counts.nx = 2;

% objective
setup.nonlin.lagrange = '2*y1^2*y2^2 + 1.25/(y2^2) + u2/y2 + u1^2 + u2^2';

% system dynamics
setup.nonlin.dynamics = '[y1 + u1/y2 + y1*y2*u2; -y2*(0.5 + y2*u2)]';

% simple bounds
UB(1).right = 4; UB(1).matrix = [1,1]; % initial states
LB(1).right = 4; LB(1).matrix = [1,1];
setup.lq.ub = UB; setup.lq.lb = LB;

% guess
Y0 = [[1,1];[1,1]];
U0 = [[0,0];[0,0]];
setup.method.guess.X = [U0,Y0];

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = HagerHouRao1_output(T,U,Y,P,F,in,opts);

% plots
HagerHouRao1_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts

% test number
num = 2;

switch num
case 1
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 50; % number of nodes
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.solver.tolerance = 1e-12;
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 7; % number of nodes
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.solver.tolerance = 1e-12;

end

end