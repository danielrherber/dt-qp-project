%--------------------------------------------------------------------------
% TwoLinkRobot.m
% Section 12.4.2 of R. Luus, Iterative Dynamic Programming. CRC Press,
% 2000, isbn: 1584881488
%--------------------------------------------------------------------------
% Similar to http://www.ee.ic.ac.uk/ICLOCS/ExampleRobotArm.html, but the
% equations at this link have some differences from Luus2000
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
y0 = [0 0 0.5 0]; % initial states
yf = [0 0 0.5 0.522]; % final states

% time horizon (scaled)
setup.t0 = 0; setup.tf = 1;

% number of controls, states, and parameters
setup.counts.nu = 2;
setup.counts.nx = 4;
setup.counts.np = 1;

% Mayer term
M(1).right = 3; M(1).left = 0; M(1).matrix = 1; % parameters
setup.lq.mayer = M;

% Lagrange term
% L(1).right = 1; L(1).left = 1; L(1).matrix = 0.01*eye(2); % controls
% setup.lq.lagrange = L;

% system dynamics
str{1} = '[';
str{2} = 'p1*( (sin(y3).*(9/4*cos(y3).*y1.^2+2*y2.^2) + 4/3*(u1-u2) - 3/2*cos(y3).*u2 )./ (31/36 + 9/4*sin(y3).^2) ); ';
str{3} = 'p1*( -( sin(y3).*(9/4*cos(y3).*y2.^2+7/2*y1.^2) - 7/3*u2 + 3/2*cos(y3).*(u1-u2) )./ (31/36 + 9/4*sin(y3).^2) ); ';
str{4} = 'p1*( y2-y1 ); ';
str{5} = 'p1*( y1 )';
str{6} = ']';
setup.nonlin.dynamics = horzcat(str{:});

% simple bounds
UB(1).right = 4; UB(1).matrix = y0; % initial states
LB(1).right = 4; LB(1).matrix = y0;
UB(2).right = 5; UB(2).matrix = yf; % final states
LB(2).right = 5; LB(2).matrix = yf;
UB(3).right = 1; UB(3).matrix = [1 1]; % controls
LB(3).right = 1; LB(3).matrix = [-1 -1];
LB(4).right = 3; LB(4).matrix = 0.1; % parameters
setup.lq.ub = UB; setup.lq.lb = LB;

% guess
Y0 = [y0;yf];
U0 = [[1 1];[-1 -1]];
P0 = [3.1;3.1];
setup.method.guess.X = [U0,Y0,P0];

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = TwoLinkRobot_output(T,U,Y,P,F,in,opts);

% plots
TwoLinkRobot_plot(T,U,Y,P,F,in,opts,sol)

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
    opts.dt.nt = 200; % number of nodes
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
end

end