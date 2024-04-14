%--------------------------------------------------------------------------
% MountainCar.m
% A. A. Melnikov, A. Makmal, H. J. Briegel, "Projective Simulation Applied
% to the Grid-World and the Mountain-Car Problem." ArXiv:1405.5459 [Cs],
% May 2014. arXiv.org, http://arxiv.org/abs/1405.5459
% Also see:
% https://openmdao.github.io/dymos/examples/mountain_car/mountain_car.html
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
setup.t0 = 0; setup.tf = 1;

% number of controls, states, and parameters
setup.counts.nu = 1;
setup.counts.nx = 2;
setup.counts.np = 1;

% Mayer term
M(1).right = 3; % parameters
M(1).left = 0; % singleton
M(1).matrix = 1;
setup.lq.mayer = M;

% system dynamics
strD{1} = '[';
strD{end+1} = 'p1*y2;';
strD{end+1} = 'p1*(0.001*u1 - 0.0025*cos(3*y1));';
strD{end+1} = ']';
setup.nonlin.dynamics = horzcat(strD{:});

% simple bounds
UB(1).right = 4; UB(1).matrix = [-0.5,0]'; % initial states
LB(1).right = 4; LB(1).matrix = [-0.5,0]';
UB(2).right = 1; UB(2).matrix = 1; % controls
LB(2).right = 1; LB(2).matrix = -1;
UB(3).right = 5; UB(3).matrix = [0.5,inf]'; % final states
LB(3).right = 5; LB(3).matrix = [0.5,0]';
UB(4).right = 3; UB(4).matrix = 10000; % parameters
LB(4).right = 3; LB(4).matrix = 0.05;
UB(5).right = 2; UB(5).matrix = [0.5,0.07]'; % states
LB(5).right = 2; LB(5).matrix = [-1.2,-0.07]';
setup.lq.ub = UB; setup.lq.lb = LB;

% guess
Y0 = [[-0.5,0];[0.5,1]];
U0 = [0;1];
P0 = [500;500];
setup.method.guess.X = [U0,Y0,P0];

% scaling
scaling(1).right = 1; % controls
scaling(1).matrix = 1;
scaling(2).right = 2; % states
scaling(2).matrix = [1.2, 0.07];
scaling(3).right = 1; % parameters
scaling(3).matrix = 100;
setup.method.scaling = scaling;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = MountainCar_output(T,U,Y,P,F,in,opts);

% plots
MountainCar_plot(T,U,Y,P,F,in,opts,sol)

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
    opts.dt.nt = 1000; % number of nodes
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.solver.maxiters = 20000;
    opts.method.derivatives = 'symbolic';

end

end