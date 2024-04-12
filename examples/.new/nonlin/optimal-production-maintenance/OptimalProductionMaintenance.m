%--------------------------------------------------------------------------
% OptimalProductionMaintenance.m
% D. I. Cho, P. L. Abad, and M. Parlar, "Optimal production and maintenance
% decisions when a system experience age-dependent deterioration," Optim.
% Control Appl. Meth., vol. 14, no. 3, pp. 153–167, Jul. 1993,
% doi: 10.1002/oca.4660140302
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
tf = 1;
a = 0;
b = 10;
h = 1;
r = 2;
q = 0;
% d = 0;
% w = 8;
c = 2.5;
alpha = 2;
s = 4;
M_ = 4;
rho = 0.1;
x0 = 3;
y0 = 1;

% time horizon
setup.t0 = 0; setup.tf = tf;

% number of controls, states, and parameters
setup.counts.nu = 2;
setup.counts.nx = 2;

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = {@(t) exp(-rho*t)*r , 0; 0, 0};
L(2).left = 0; L(2).right = 1; L(2).matrix = {@(t) exp(-rho*t)*q, @(t) exp(-rho*t)*c};
L(3).left = 0; L(3).right = 2; L(3).matrix = {@(t) exp(-rho*t)*h, 0};
setup.lq.lagrange = L;

% Mayer term
M(1).left = 0; M(1).right = 5; M(1).matrix = -[a*exp(-rho*tf), b*exp(-rho*tf)]; % final states
setup.lq.mayer = M;

% system dynamics
str = {};
str{end+1} = '[';
str{end+1} = 'y2*u1 - s;';
str{end+1} = '-(alpha_ + u2)*y2 + u2;';
str{end+1} = ']';
str = horzcat(str{:});
setup.nonlin.dynamics = str;

% symbolic data for nonlin
setup.nonlin.data.symbols = 'alpha_ s';
setup.nonlin.data.values = [alpha s];

% simple bounds
UB(1).right = 4; UB(1).matrix = [x0,y0]; % initial states
LB(1).right = 4; LB(1).matrix = [x0,y0];
UB(2).right = 1; UB(2).matrix = [inf,M_]; % controls
LB(2).right = 1; LB(2).matrix = [0,0];
LB(3).right = 2; LB(3).matrix = [0,-inf];
setup.lq.ub = UB; setup.lq.lb = LB;

% guess
Y0 = [[x0,y0];[0,y0]];
U0 = [[0,0];[0,M_]];
setup.method.guess.X = [U0,Y0];

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = OptimalProductionMaintenance_output(T,U,Y,P,F,in,opts);

% plots
OptimalProductionMaintenance_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts

% test number
num = 1;

switch num
    case 1
    % default parameters
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 300; % number of nodes
    opts.solver.function = 'ipfmincon';
    opts.solver.tolerance = 1e-8;
    opts.method.form = 'nonlinearprogram';
end

end