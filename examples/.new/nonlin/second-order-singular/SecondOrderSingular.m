%--------------------------------------------------------------------------
% SecondOrderSingular.m
% G. M. Aly, "The computation of optimal singular control", International
% Journal of Control, vol. 28, no. 5, pp. 681–688, Nov. 1978,
% doi: 10.1080/00207177808922489
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
setup.t0 = 0; setup.tf = 5;

% number of controls, states, and parameters
setup.counts.nu = 1;
setup.counts.nx = 3;

% Mayer term
M(1).left = 0; M(1).right = 5; M(1).matrix = [0 0 1]; % final states
setup.lq.mayer = M;

% system dynamics
str{1} = '[';
str{end+1} = 'y2; ';
str{end+1} = 'u1; ';
str{end+1} = 'y1^2/2 + y2^2/2';
str{end+1} = ']';
setup.nonlin.dynamics = horzcat(str{:});

% simple bounds
UB(1).right = 4; UB(1).matrix = [0,1,0]; % initial states
LB(1).right = 4; LB(1).matrix = [0,1,0];
UB(2).right = 1; UB(2).matrix = 1; % controls
LB(2).right = 1; LB(2).matrix = -1;
setup.lq.ub = UB; setup.lq.lb = LB;

% guess
Y0 = [[0,1,0];[0,1,1]];
U0 = [-1;-1];
setup.method.guess.X = [U0,Y0];

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = SecondOrderSingular_output(T,U,Y,P,F,in,opts);

% plots
SecondOrderSingular_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts

% test number
num = 1;

switch num
case 1
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100; % number of nodes
    opts.solver.maxiters = 300;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.method.olqflag = true;
    opts.method.derivatives = 'complex';
    opts.solver.tolerance = 1e-12;
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 100; % number of nodes
    opts.solver.maxiters = 300;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.method.olqflag = true;
    opts.method.derivatives = 'complex';
end

end