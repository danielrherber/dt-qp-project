%--------------------------------------------------------------------------
% Vanderpol.m
% E. B. Canto, et al., "Restricted second order information for the
% solution of optimal control problems using control vector
% parameterization", *Journal of Process Control*, vol. 12, no. 2002,
% pp. 243-255, 2002, doi: 10.1016/S0959-1524(01)00008-7
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
setup = createSetup(opts);

% solve with DTQP
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

% post-processing
[O,sol] = postProcessing(T,U,Y,P,F,in,opts);

%--------------------------------------------------------------------------
% create setup structure
%--------------------------------------------------------------------------
function setup = createSetup(opts)

% initialize struct
setup = DTQP_setup_initialize;

% time horizon
setup.t0 = 0; setup.tf = 5;

% select case (see below)
casenum = 2;

switch casenum
    %----------------------------------------------------------------------
    case 1

    % number of controls, states, and parameters
    setup.counts.nu = 1; % number of controls
    setup.counts.nx = 2; % number of states

    % Lagrange term
    setup.nonlin.lagrange = 'y1^2 + y2^2 + u1^2';
    % L(1).left = 1; L(1).right = 1; L(1).matrix = 1;
    % L(2).left = 2; L(2).right = 2; L(2).matrix = eye(2);
    % setup.lq.lagrange = lq;

    % system dynamics
    setup.nonlin.dynamics = '[y2;-y1+y2-y1^2*y2+u1]';

    % simple bounds
    UB(1).right = 4; UB(1).matrix = [1;0]; % initial states
    LB(1).right = 4; LB(1).matrix = [1;0];
    UB(2).right = 1; UB(2).matrix = 1; % controls
    LB(2).right = 1; LB(2).matrix = -0.3;
    setup.lq.ub = UB; setup.lq.lb = LB;

    % guess
    Y0 = [[1,0];[1,0]];
    U0 = [0;0];
    setup.method.guess.X = [U0,Y0];

    %----------------------------------------------------------------------
    case 2 % control co-design example

    % auxiliary data
    dmin = [0.1 0.1];
    dmax = [5 5];

    % number of controls, states, and parameters
    setup.counts.nu = 1; % number of controls
    setup.counts.nx = 2; % number of states
    setup.counts.np = 2; % number of parameters

    % Lagrange term
    if ~isfield(opts.method,'olqflag') || opts.method.olqflag
        L(1).left = 1; L(1).right = 1; L(1).matrix = 1;
        L(2).left = 2; L(2).right = 2; L(2).matrix = eye(2);
        setup.lq.lagrange = L;
    else
        setup.nonlin.lagrange = 'y1^2 + y2^2 + u1^2';
    end

    % system dynamics
    setup.nonlin.dynamics = '[y2;-p1*y1+y2-p2*y1^2*y2+u1]';

    % simple bounds
    UB(1).right = 4; UB(1).matrix = [1;0]; % initial states
    LB(1).right = 4; LB(1).matrix = [1;0]; % initial states
    UB(2).right = 1; UB(2).matrix = 1; % controls
    LB(2).right = 1; LB(2).matrix = -0.5; % controls
    UB(3).right = 3; UB(3).matrix = dmax; % parameters
    LB(3).right = 3; LB(3).matrix = dmin; % parameters
    LB(4).right = 2; LB(4).matrix = [-inf;-0.4]; % states
    setup.lq.ub = UB; setup.lq.lb = LB;

    % guess
    Y0 = [[1,0];[0,0]];
    U0 = [0;0];
    P0 = [[1,1];[1,1]];
    setup.method.guess.X = [U0,Y0,P0];

    %----------------------------------------------------------------------
    case 3 % completely symbolic formulation

    % number of controls, states, and parameters
    setup.counts.nu = 1; % number of controls
    setup.counts.nx = 2; % number of states

    % Lagrange term
    setup.nonlin.lagrange = 'y1^2 + y2^2 + u1^2';

    % system dynamics
    setup.nonlin.dynamics = '[y2; -y1 + y2 - y1^2*y2 + u1]';

    % equality constraints
    setup.nonlin.equality.func = '[yi1 - 1; yi2 - 0]';

    % inequality constraints
    setup.nonlin.inequality.func = '[u1 - 1; -u1 - 0.3]';

    % guess
    Y0 = [[1,0];[1,0]];
    U0 = [0;0];
    setup.method.guess.X = [U0,Y0];

end

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = Vanderpol_output(T,U,Y,P,F,in,opts);

% plots
Vanderpol_plot(T,U,Y,P,F,in,opts,sol)

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
    opts.method.form = 'nonlinearprogram';
    opts.solver.function = 'ipfmincon';
    opts.solver.maxiters = 4000;
    opts.solver.display = 'iter';
    opts.method.olqflag = true;
    opts.method.derivativeflag = true;
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 40; % number of nodes
    opts.method.form = 'nonlinearprogram';
    opts.solver.function = 'ipfmincon';
    opts.solver.display = 'iter';
    opts.method.olqflag = true;
    opts.method.derivativeflag = true;
case 3
    opts.general.displevel = 1;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.method.maxiters = 5000;
    opts.dt.nt = 40; % number of nodes
    opts.method.form = 'qlin';
    opts.solver.display = 'none';
    opts.method.trustregionflag = false;
    opts.method.improveguess = false;
end

end