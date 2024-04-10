%--------------------------------------------------------------------------
% FreeFlyingRobot.m
% 326-330 of J. T. Betts, "Practical Methods for Optimal Control and
% Estimation Using Nonlinear Programming*." Society for Industrial and
% Applied Mathematics, Jan. 2010, doi: 10.1137/1.9780898718577
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
tf = 12;
a = 0.2; b = 0.2;

% time horizon
setup.t0 = 0; setup.tf = tf;

% number of controls, states, and parameters
setup.counts.nu = 4;
setup.counts.nx = 6;

% Lagrange term
% setup.nonlin.lagrange = 'u1+u2+u3+u4';
L(1).left = 0; L(1).right = 1; L(1).matrix = [1 1 1 1];
setup.lq.lagrange = L;

% system dynamics
str{1} = '[';
str{end+1} = 'y4;';
str{end+1} = 'y5;';
str{end+1} = 'y6;';
str{end+1} = '(u1-u2+u3-u4)*cos(y3);';
str{end+1} = '(u1-u2+u3-u4)*sin(y3);';
str{end+1} = 'a*(u1-u2) - b*(u3-u4)';
str{end+1} = ']';
setup.nonlin.dynamics = horzcat(str{:});

% data for symbolic functions
setup.nonlin.data.symbols = 'a b';
setup.nonlin.data.values = [a b];

% inequality constraints
inequality.func = '[u1+u2-1;u3+u4-1]';
inequality.pathboundary = [1 1];
setup.nonlin.inequality = inequality;

% simple bounds
UB(1).right = 4; UB(1).matrix = [-10;-10;pi/2;0;0;0]; % initial states
LB(1).right = 4; LB(1).matrix = [-10;-10;pi/2;0;0;0];
UB(2).right = 5; UB(2).matrix = [0;0;0;0;0;0]; % final states
LB(2).right = 5; LB(2).matrix = [0;0;0;0;0;0];
UB(3).right = 1; UB(3).matrix = 0.3+[1;1;1;1]; % controls
LB(3).right = 1; LB(3).matrix = [0;0;0;0];
setup.lq.ub = UB; setup.lq.lb = LB;

% guess
Y0 = [[-10,-10,pi/2,0,0,0];[0,0,0,0,0,0]];
U0 = [[0,0,0,0];[0,0,0,0]];
setup.method.guess.X = [U0,Y0];

% scaling
scaling(1).right = 1; % controls
scaling(1).matrix = [1;1;1;1];
scaling(2).right = 2; % states
scaling(2).matrix = [10,10,pi/2,1,1,1];
setup.scaling = scaling;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = FreeFlyingRobot_output(T,U,Y,P,F,in,opts);

% plots
FreeFlyingRobot_plot(T,U,Y,P,F,in,opts,sol)

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
    opts.dt.nt = 300; % number of nodes
    opts.solver.display = 'iter';
    opts.solver.maxiters = 2000;
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.method.olqflag = true;
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 100; % number of nodes
    opts.solver.display = 'iter';
    opts.solver.maxiters = 2000;
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.method.olqflag = true;
end

end