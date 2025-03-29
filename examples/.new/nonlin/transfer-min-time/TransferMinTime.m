%--------------------------------------------------------------------------
% TransferMinTime.m
% p. 66-69 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
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
xT  = 0.1405;
m0 = 1;
dm = 0.1405*0.533;
xmu = 1;
rf = 1.525;

% time horizon
setup.t0 = 0; setup.tf = 1;

% number of controls, states, and parameters
setup.counts.nu = 1;
setup.counts.nx = 4;
setup.counts.np = 1;

% Mayer term
M(1).left = 0; M(1).right = 3; M(1).matrix = 1;
setup.lq.mayer = M;

% system dynamics
str{1} = '[';
str{end+1} = 'p1*(y3); ';
str{end+1} = 'p1*(y4/y1); ';
str{end+1} = 'p1*(y4^2/y1 - xmu/y1^2 + (xT/(m0 - dm*p1*t))*sin(u1)); ';
str{end+1} = 'p1*(-y3*y4/y1 + (xT/(m0 - dm*p1*t))*cos(u1))';
str{end+1} = ']';
setup.nonlin.dynamics = horzcat(str{:});

% symbolic data for nonlin
setup.nonlin.data.symbols = 'xmu xT m0 dm';
setup.nonlin.data.values = [xmu xT m0 dm];

% simple bounds
UB(1).right = 4; UB(1).matrix = [1,0,0,1]; % initial states
LB(1).right = 4; LB(1).matrix = [1,0,0,1];
UB(2).right = 5; UB(2).matrix = [rf,inf,0,sqrt(xmu/rf)]; % final states
LB(2).right = 5; LB(2).matrix = [rf,-inf,0,sqrt(xmu/rf)];
UB(3).right = 3; UB(3).matrix = [100]; % parameters
LB(3).right = 3; LB(3).matrix = [0];
UB(4).right = 1; UB(4).matrix = [2*pi]; % controls
LB(4).right = 1; LB(4).matrix = [0];
setup.lq.ub = UB; setup.lq.lb = LB;

% guess
Y0 = [[1,0,0,1];[rf,pi,0,sqrt(xmu/rf)]];
U0 = [[0];[2*pi]];
P0 = [[3];[3]];
setup.method.guess.X = [U0,Y0,P0];

% scaling
scaling(1).right = 1; % controls
scaling(1).matrix = [pi];
scaling(2).right = 2; % states
scaling(2).matrix = [1.5,2.5,0.3,1];
setup.method.scaling = scaling;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = TransferMinTime_output(T,U,Y,P,F,in,opts);

% plots
TransferMinTime_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts

% test number
num = 2;

switch num
case 1 % ipfmincon method
    opts.general.plotflag = 1; % create the plots
    opts.general.displevel = 2;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 40; % number of nodes
    opts.solver.tolerance = 1e-8;
    opts.solver.maxiters = 2000;
    opts.solver.display = 'iter';
    opts.method.form = 'nonlinearprogram';
case 2
    opts.general.plotflag = 1; % create the plots
    opts.general.displevel = 2;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 20; % number of nodes
    opts.solver.tolerance = 1e-8;
    opts.solver.maxiters = 2000;
    opts.solver.display = 'iter';
    opts.method.form = 'nonlinearprogram';
    opts.dt.meshr.method = 'RICHARDSON-DOUBLING';
    opts.dt.meshr.tolerance = 1e-6;
case 3
    opts.general.plotflag = 1; % create the plots
    opts.general.displevel = 2;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 10; % number of nodes
    opts.solver.tolerance = 1e-8;
    opts.solver.display = 'iter';
    opts.method.form = 'nonlinearprogram';
    opts.dt.meshr.method = 'RICHARDSON-DOUBLING';
    opts.dt.meshr.tolerance = 1e-6;
end

end