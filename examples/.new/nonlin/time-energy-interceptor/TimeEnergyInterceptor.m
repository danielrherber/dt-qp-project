%--------------------------------------------------------------------------
% TimeEnergyInterceptor.m
% A. Banerjee, M. Nabi, and T. Raghunathan, "Time-energy optimal guidance
% strategy for realistic interceptor using pseudospectral method",
% Transactions of the Institute of Measurement and Control, vol. 42,
% no. 13. SAGE Publications, pp. 2361-2371, Apr. 21, 2020
% doi:10.1177/0142331220910919
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
engagement_case = 1; % (see below)

switch engagement_case

case 1 % tail chase (nominal)
    VT = 200; % target speed [m/s]
    VM = 1000; % missile speed [m/s]
    aT = @(t) zeros(size(t)); % lateral acceleration of the target [m/s^2]
    AT0 = 0; % initial target flight path angle [rad]
    RT10 = 20000; % initial target position along axis 1 [m/s]
    RT20 = 10000; % initial target position along axis 1 [m/s]
    AM0 = 17*pi/18; % initial missile flight path angle [rad]
    RM10 = 0; % initial missile position along axis 1 [m/s]
    RM20 = 10000; % initial missile position along axis 2 [m/s]
    Wt = 10000; % time weight in the objective function
    ambar = 150; % control input limit [m/s^2]
    r = 5; % lethal radius [m]

case 3 % head-on (nominal)
    VT = 200; % target speed [m/s]
    VM = 1000; % missile speed [m/s]
    aT = @(t) zeros(size(t)); % lateral acceleration of the target [m/s^2]
    AT0 = pi; % initial target flight path angle [rad]
    RT10 = 20000; % initial target position along axis 1 [m/s]
    RT20 = 10000; % initial target position along axis 1 [m/s]
    AM0 = 17*pi/18; % initial missile flight path angle [rad]
    RM10 = 0; % initial missile position along axis 1 [m/s]
    RM20 = 10000; % initial missile position along axis 2 [m/s]
    Wt = 10000; % time weight in the objective function
    ambar = 165; % control input limit [m/s^2]
    r = 5; % lethal radius [m]
end
auxdata.r = r; % save for later

% time horizon
setup.t0 = 0; setup.tf = 1;

% symbolic data for nonlin
setup.nonlin.data.symbols = 'VT VM aT Wt r';
setup.nonlin.data.values = {VT VM aT Wt r};

% Lagrange term
setup.nonlin.lagrange = 'p1*(u1^2/Wt + 1)';

% select the formulation to use
formulation = 2;
auxdata.formulation = formulation; % save for later

switch formulation
    %----------------------------------------------------------------------
    case 1 % original

    % number of controls, states, and parameters
    setup.counts.nu = 1;
    setup.counts.nx = 6;
    setup.counts.np = 1;

    % system dynamics
    % u1 = am, y1 = AT, y2 = RT1, y3 = RT2, y4 = AM, y5 = RM1, y6 = RM2, p1 = tf
    strD{1} = '[';
    strD{end+1} = 'p1*(aT/VT);';
    strD{end+1} = 'p1*(VT*cos(y1));';
    strD{end+1} = 'p1*(VT*sin(y1));';
    strD{end+1} = 'p1*(u1/VM);';
    strD{end+1} = 'p1*(VM*cos(y4));';
    strD{end+1} = 'p1*(VM*sin(y4));';
    strD{end+1} = ']';
    setup.nonlin.dynamics = horzcat(strD{:});

    % equality constraints
    setup.nonlin.equality.func = '[(yf2-yf5)^2/r^2 + (yf3-yf6)^2/r^2 - 1]';
    setup.nonlin.equality.pathboundary = 0;

    % initial states
    Y0_ = [AT0 RT10 RT20 AM0 RM10 RM20];

    % guess
    U0 = [ambar;-ambar];
    P0 = [50;50];
    Y0 = [Y0_;[AT0 RT10+P0(1)^2/2 RT20+0 0 0 0]];
    setup.method.guess.X = [U0,Y0,P0];

    % scaling
    scaling(1).right = 1; % controls
    scaling(1).matrix = ambar;
    scaling(2).right = 2; % states
    scaling(2).matrix = [pi 1e4 1e4 pi 1e4 1e4];
    scaling(3).right = 3; % parameters
    scaling(3).matrix = 100;
    setup.method.scaling = scaling;

    %----------------------------------------------------------------------
    case 2 % reduced formulation

    % number of controls, states, and parameters
    setup.counts.nu = 1;
    setup.counts.nx = 4;
    setup.counts.np = 1;

    % system dynamics
    % u1 = am, y1 = AT, y2 = RT1, y3 = RT2, y4 = AM, y5 = RM1, y6 = RM2, p1 = tf
    strD{1} = '[';
    strD{end+1} = 'p1*(aT/VT);';
    strD{end+1} = 'p1*(u1/VM);';
    strD{end+1} = 'p1*(VT*cos(y1)-VM*cos(y2));';
    strD{end+1} = 'p1*(VT*sin(y1)-VM*sin(y2));';
    strD{end+1} = ']';
    setup.nonlin.dynamics = horzcat(strD{:});

    % Lagrange term
    setup.nonlin.lagrange = 'p1*(u1^2/Wt + 1)';

    % equality constraints
    setup.nonlin.equality.func = '[yf3^2/r^2 + yf4^2/r^2 - 1]';
    setup.nonlin.equality.pathboundary = 0;

    % initial states
    Y0_ = [AT0 AM0 RT10-RM10 RT20-RM20];

    % guess
    U0 = [ambar;-ambar];
    P0 = [50;50];
    Y0 = [Y0_;[0 0 RT10+P0(1)^2/2 RT20]];
    setup.method.guess.X = [U0,Y0,P0];

    % scaling
    scaling(1).right = 1; % controls
    scaling(1).matrix = ambar;
    scaling(2).right = 2; % states
    scaling(2).matrix = [pi pi 1e4 1e4];
    scaling(3).right = 3; % parameters
    scaling(3).matrix = 100;
    setup.method.scaling = scaling;

end

% simple bounds
UB(1).right = 4; UB(1).matrix = Y0_; % initial states
LB(1).right = 4; LB(1).matrix = Y0_;
UB(2).right = 1; UB(2).matrix = ambar; % controls
LB(2).right = 1; LB(2).matrix = -ambar;
UB(3).right = 3; UB(3).matrix = 200; % parameters
LB(3).right = 3; LB(3).matrix = 10;
setup.lq.ub = UB; setup.lq.lb = LB;

% assign
setup.auxdata = auxdata;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = TimeEnergyInterceptor_output(T,U,Y,P,F,in,opts);

% plots
TimeEnergyInterceptor_plot(T,U,Y,P,F,in,opts,sol)

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
    opts.dt.nt = 10000; % number of nodes
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.method.derivatives = 'symbolic';
    opts.solver.tolerance = 1e-10;
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 40; % number of nodes
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.method.derivatives = 'symbolic';
    opts.solver.tolerance = 1e-10;
end

end