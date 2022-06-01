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
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = TimeEnergyInterceptor(varargin)
% input arguments can be provided in the format 'TimeEnergyInterceptor(auxdata,opts)'

% set local functions
ex_opts = @TimeEnergyInterceptor_opts; % options function
ex_output = @TimeEnergyInterceptor_output; % output function
ex_plot = @TimeEnergyInterceptor_plot; % plot function

% set auxdata and opts (see local_opts)
[auxdata,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
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

%% setup
% time horizon
auxdata.t0 = 0; auxdata.tf = 1;

% number of controls, states, and parameters
n.nu = 1; n.ny = 6; n.np = 1;

% problem parameters
element.parameter_list = 'VT VM aT Wt r';
element.parameter_values = {VT VM aT Wt r};

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
element.dynamics = horzcat(strD{:});

% Lagrange term
element.lagrange = 'p1*(u1^2/Wt + 1)';

% equality constraints
element.h.func = '[(yf2-yf5)^2/r^2 + (yf3-yf6)^2/r^2 - 1]';
element.h.pathboundary = 0;

% simple bounds
Y0_ = [AT0 RT10 RT20 AM0 RM10 RM20];
UB(1).right = 4; UB(1).matrix = Y0_; % initial states
LB(1).right = 4; LB(1).matrix = Y0_;
UB(2).right = 1; UB(2).matrix = ambar; % controls
LB(2).right = 1; LB(2).matrix = -ambar;
UB(3).right = 3; UB(3).matrix = 200; % parameters
LB(3).right = 3; LB(3).matrix = 10;

% guess
U0 = [[0];[0]];
P0 = [[100];[100]];
Y0 = [Y0_;[AT0 RT10+P0(1)^2/2 RT20+0 0 0 0]];
setup.guess.X = [U0,Y0,P0];

% scaling
setup.scaling(1).right = 1; % controls
setup.scaling(1).matrix = ambar;
setup.scaling(2).right = 2; % states
setup.scaling(2).matrix = [pi 1e4 1e4 pi 1e4 1e4];
setup.scaling(3).right = 3; % parameters
setup.scaling(3).matrix = 100;

% combine structures
setup.element = element; setup.UB = UB; setup.LB = LB;
setup.t0 = auxdata.t0; setup.tf = auxdata.tf; setup.auxdata = auxdata; setup.n = n;

%% solve
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = ex_output(T,U,Y,P,F,in,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
ex_plot(T,U,Y,P,F,in,opts,sol)

end
% User options function for this example
function opts = TimeEnergyInterceptor_opts
% test number
num = 2;

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