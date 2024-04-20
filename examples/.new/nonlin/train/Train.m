%--------------------------------------------------------------------------
% Train.m
% R. J. Vanderbei, "Case Studies in Trajectory Optimization: Trains,
% Planes, and Other Pastimes," Optimization and Engineering, vol. 2, no. 2,
% pp. 215-243, 2001, doi: 10.1023/a:1013145328012
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
tf = 4.8;
a = 0.3; b = 0.14; c = 0.16;
z1 = 2; z2 = 4;
s1 = 2; s2 = 0; s3 = -2;

% time horizon
setup.t0 = 0; setup.tf = tf;

% number of controls, states, and parameters
setup.counts.nu = 2;
setup.counts.nx = 2;

% Lagrange term
% setup.nonlin.lagrange = 'u1*y2 + 10^-3*(u1^2 + u2^2)';
L(1).left = 1; L(1).right = 1; L(1).matrix = diag([1e-3,1e-3]);
L(2).left = 1; L(2).right = 2; L(2).matrix = [0,0;1,0;0,0;0,0];
setup.lq.lagrange = L;

% system dynamics
str{1} = '[';
str{end+1} = 'y2;';
str{end+1} = '(s2-s1)/pi*atan((y1-2)/0.05)+(s3-s2)/pi*atan((y1-4)/0.05)-(a+b*y2+c*y2^2)+u1-u2';
str{end+1} = ']';
setup.nonlin.dynamics = horzcat(str{:});

% symbolic data for nonlin
setup.nonlin.data.symbols = 'a b c z1 z2 s1 s2 s3';
setup.nonlin.data.values = [a b c z1 z2 s1 s2 s3];

% simple bounds
UB(1).right = 4; UB(1).matrix = [0;0]; % initial states
LB(1).right = 4; LB(1).matrix = [0;0];
UB(2).right = 5; UB(2).matrix = [6;0]; % final states
LB(2).right = 5; LB(2).matrix = [6;0];
UB(3).right = 1; UB(3).matrix = [10;2]; % controls
LB(3).right = 1; LB(3).matrix = [0;0];
setup.lq.ub = UB; setup.lq.lb = LB;

% guess
Y0 = [[0,0];[6,0]];
U0 = [[10,2];[0,0]];
setup.method.guess.X = [U0,Y0];

% scaling
scaling(1).right = 1; % controls
scaling(1).matrix = [10,2];
scaling(2).right = 2; % states
scaling(2).matrix = [6,6];
setup.method.scaling = scaling;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = Train_output(T,U,Y,P,F,in,opts);

% plots
Train_plot(T,U,Y,P,F,in,opts,sol)

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
    opts.dt.nt = 500; % number of nodes
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.method.olqflag = true;
end

end