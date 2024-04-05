%--------------------------------------------------------------------------
% Brachistochrone.m
% pp. 133-134 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
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

% select problem variant (see below)
auxdata.variant = 4; 

switch auxdata.variant
    case 1 % final x and y specified
    auxdata.g = 10;
    auxdata.xf = 2;
    auxdata.yf = 2;
    case 2 % final x and y specified
    auxdata.g = 10;
    auxdata.xf = 10;
    auxdata.yf = 2;
    case 3 % final x specified
    auxdata.g = 32.1740;
    auxdata.xf = 1;
    case 4 % final x specified with path constraint
    auxdata.g = 32.1740;
    auxdata.xf = 1;
    auxdata.h = 0.1;
    auxdata.theta = atan(0.5);
end

% other auxiliary data
auxdata.t0 = 0; auxdata.tf = 1;
setup.auxdata = auxdata;

% time horizon
setup.t0 = auxdata.t0; setup.tf = auxdata.tf;

% number of controls, states, and parameters
counts.nu = 1;
counts.nx = 3;
counts.np = 1;
counts.nv = 0;
setup.counts = counts;

% Mayer term
M(1).left = 0; M(1).right = 3; M(1).matrix = 1; % parameters
setup.lq.mayer = M;

% simple bounds
UB(1).right = 4; UB(1).matrix = [0,0,0]; % initial states
LB(1).right = 4; LB(1).matrix = [0,0,0];
UB(3).right = 1; UB(3).matrix = pi; % controls
LB(3).right = 1; LB(3).matrix = -pi;
UB(4).right = 3; UB(4).matrix = 200; % parameters
LB(4).right = 3; LB(4).matrix = 0;

switch auxdata.variant
    %----------------------------------------------------------------------
    case {1,2}

    % system dynamics
    setup.nonlin.dynamics = '[p1*y3*sin(u1); p1*y3*cos(u1); p1*g*cos(u1)]';
    setup.nonlin.data.symbols = 'g';
    setup.nonlin.data.values = [auxdata.g];

    % simple bounds
    UB(2).right = 5; UB(2).matrix = [auxdata.xf,auxdata.yf,inf]; % final states
    LB(2).right = 5; LB(2).matrix = [auxdata.xf,auxdata.yf,-inf];
    setup.lq.ub = UB; setup.lq.lb = LB;

    % guess
    Y0 = [[0,0,0];[auxdata.xf,auxdata.yf,0]];
    U0 = [0;0];
    P0 = [1;1];
    setup.method.guess.X = [U0,Y0,P0];

    %----------------------------------------------------------------------
    case {3,4}

    % system dynamics
    setup.nonlin.dynamics = '[p1*y3*cos(u1); p1*y3*sin(u1); p1*g*sin(u1)]';
    setup.nonlin.data.symbols = 'g';
    setup.nonlin.data.values = [auxdata.g];

    % simple bounds
    UB(2).right = 5; UB(2).matrix = [auxdata.xf,inf,inf]; % final states
    LB(2).right = 5; LB(2).matrix = [auxdata.xf,-inf,-inf];
    setup.lq.ub = UB; setup.lq.lb = LB;

    % guess
    Y0 = [[0,0,0];[auxdata.xf,0,0]];
    U0 = [0;0];
    P0 = [1;1];
    setup.method.guess.X = [U0,Y0,P0];

end

% state inequality constraint
if auxdata.variant == 4
    inequality.func = 'y2 - y1/2 - 0.1'; % hard-coded at the moment
    inequality.pathboundary = 1;
    setup.nonlin.inequality = inequality;
end

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = Brachistochrone_output(T,U,Y,P,F,in,opts);

% plots
Brachistochrone_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts

% test number
num = 1;

switch num
case 1 % ipfmincon method
    opts.general.plotflag = 1; % create the plots
    opts.general.saveflag = false;
    opts.general.displevel = 2;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 50; % number of nodes
    opts.solver.tolerance = 1e-11;
    opts.solver.maxiters = 1000;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 50; % number of nodes
    opts.solver.tolerance = 1e-11;
    opts.solver.maxiters = 100;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
end

end