%--------------------------------------------------------------------------
% AlpRider.m
% pp. 163-165 in J. T. Betts, "Practical Methods for Optimal Control and
% Estimation Using Nonlinear Programming." Society for Industrial and
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
auxdata.t0 = 0; auxdata.tf = 20;
setup.auxdata = auxdata;

% time horizon
setup.t0 = auxdata.t0; setup.tf = auxdata.tf;

% number of controls, states, and parameters
setup.counts.nu = 2;
setup.counts.nx = 4;
setup.counts.np = 0;

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = diag([1e-2,1e-2]);
L(2).left = 2; L(2).right = 2; L(2).matrix = diag([1e2,1e2,1e2,1e2]);
setup.lq.lagrange = L;

% system dynamics
str{1} = '[';
str{end+1} = '-10*y1 + u1 + u2; ';
str{end+1} = '-2*y2 + u1 + 2*u2; ';
str{end+1} = '-3*y3 + 5*y4 + u1 - u2; ';
str{end+1} = '5*y3 - 3*y4 + u1 + 3*u2';
str{end+1} = ']';
setup.nonlin.dynamics = horzcat(str{:});

% nonlinear inequality constraint
inequality.func = '3*exp(-12*(t-3)^2) + 3*exp(-10*(t-6)^2) + 3*exp(-6*(t-10)^2) + 8*exp(-4*(t-15)^2) + 0.01 - y1^2 - y2^2 - y3^2 - y4^2';
inequality.pathboundary = 1;
setup.nonlin.inequality = inequality;

% simple bounds
UB(1).right = 4; UB(1).matrix = [2,1,2,1]; % initial states
LB(1).right = 4; LB(1).matrix = [2,1,2,1];
UB(2).right = 5; UB(2).matrix = [2,3,1,-2]; % final states
LB(2).right = 5; LB(2).matrix = [2,3,1,-2];
setup.lq.ub = UB; setup.lq.lb = LB;

% guess
Y0 = [[2,1,2,1];[2,3,1,-2]];
U0 = [[0,0];[0,0]];
setup.method.guess.X = [U0,Y0];

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = AlpRider_output(T,U,Y,P,F,in,opts);

% plots
AlpRider_plot(T,U,Y,P,F,in,opts,sol)

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
        opts.dt.nt = 100; % number of nodes
        opts.solver.maxiters = 300;
        opts.solver.display = 'iter';
        opts.solver.function = 'ipfmincon';
        opts.method.form = 'nonlinearprogram';
        opts.method.olqflag = false;
        opts.dt.meshr.method = 'RICHARDSON-DOUBLING';
        opts.dt.meshr.tolerance = 1e-3;
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
        opts.method.olqflag = false;
end

end