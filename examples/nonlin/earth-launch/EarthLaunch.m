%--------------------------------------------------------------------------
% EarthLaunch.m
% pp. 222-226 in J. M. Longuski, J. J. Guzmán, and J. E. Prussing, Optimal
% Control with Aerospace Applications. Springer New York, 2014 [Online].
% Available: http://dx.doi.org/10.1007/978-1-4614-8945-0
% Values based on the following implementation of this problem:
% https://openmdao.github.io/dymos/examples/ssto_earth/ssto_earth.html
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = EarthLaunch(varargin)
% input arguments can be provided in the format 'EarthLaunch(auxdata,opts)'

% set local functions
ex_opts = @EarthLaunch_opts; % options function
ex_output = @EarthLaunch_output; % output function
ex_plot = @EarthLaunch_plot; % plot function

% set auxdata and opts (see local_opts)
[auxdata,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
% initial conditions
x0 = 0; % initial range [m]
y0 = 0; % initial altitude [m]
vx0 = 0; % initial x velocity [m/s]
vy0 = 0; % initial y velocity [m/s]
m0 = 117000; % initial mass [kg]

% final conditions
yf = 1.85e5; % final altitude [m]
vxf = 7796.6961; % final x velocity [m/s]
vyf = 0; % final y velocity [m/s]

% other parameters
g = 9.80665; % gravitational acceleration [m/s^]
rho_ref = 1.225; % reference atmospheric density [kg/m^3]
h_ref = 8.44e3; % reference altitude [m]
CD = 0.5; % coefficient of drag []
S = 7.069; % aerodynamic reference area [m^2]
FT = 2100000; % thrust [N]
Isp = 265.2; % specific impulse [s]
hs = 1e5; % distance scale factor [m]
vs = 1e3; % velocity scale factor [m/s]
ts = 100; % time scale factor [s]

%% setup
% time horizon
auxdata.t0 = 0; auxdata.tf = 1;

% number of controls, states, and parameters
n.nu = 1; n.ny = 5; n.np = 1;

% system dynamics
strD{1} = '[';
strD{end+1} = 'p1*(y3);'; % x, range [m]
strD{end+1} = 'p1*(y4);'; % y, altitude [m]
strD{end+1} = 'p1*((FT*cos(u1) - 0.5*CDA*(rho_ref*exp(-y2/hs))*y3^2)/y5);'; % vx, x velocity [m/s]
strD{end+1} = 'p1*((FT*sin(u1) - 0.5*CDA*(rho_ref*exp(-y2/hs))*y4^2)/y5 - g);'; % vy, y velocity [m/s]
strD{end+1} = 'p1*(-FT/(g*Isp));'; % m, mass [kg]
strD{end+1} = ']';
element.dynamics = horzcat(strD{:});
element.parameter_list = 'FT CDA g Isp rho_ref hs';
element.parameter_values = [FT CD*S g Isp rho_ref h_ref];

% Mayer term
M(1).right = 3; % parameters
M(1).left = 0; % singleton
M(1).matrix = 1;

% simple bounds
UB(1).right = 4; UB(1).matrix = [x0,y0,vx0,vy0,m0]'; % initial states
LB(1).right = 4; LB(1).matrix = [x0,y0,vx0,vy0,m0]';
UB(2).right = 1; UB(2).matrix = pi/2; % controls
LB(2).right = 1; LB(2).matrix = -pi/2;
UB(3).right = 5; UB(3).matrix = [inf,yf,vxf,vyf,inf]'; % final states
LB(3).right = 5; LB(3).matrix = [-inf,yf,vxf,vyf,-inf]';
UB(4).right = 3; UB(4).matrix = 500; % parameters
LB(4).right = 3; LB(4).matrix = 10;

% guess
Y0 = [[x0,y0,vx0,vy0,m0];[1.15e5,yf,vxf,vyf,1163]];
U0 = [1.5;-0.76];
P0 = [150;150];
setup.guess.X = [U0,Y0,P0];

% scaling
setup.scaling(1).right = 1; % controls
setup.scaling(1).matrix = pi/2;
setup.scaling(2).right = 2; % states
setup.scaling(2).matrix = [hs,hs,vs,vs,m0];
setup.scaling(3).right = 1; % parameters
setup.scaling(3).matrix = ts;

% combine structures
setup.element = element; setup.M = M; setup.UB = UB; setup.LB = LB;
setup.t0 = auxdata.t0; setup.tf = auxdata.tf; setup.auxdata = auxdata; setup.n = n;

%% solve
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = ex_output(T,U,Y,P,F,in,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
% disp("paused"); pause % for quasilinearization plots
ex_plot(T,U,Y,P,F,in,opts,sol)

end
% User options function for this example
function opts = EarthLaunch_opts
% test number
num = 3;

switch num
    case 1
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 400; % number of nodes
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.solver.tolerance = 1e-10;
    opts.method.form = 'nonlinearprogram';
    opts.method.derivatives = 'symbolic';

    case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 100; % number of nodes
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.solver.tolerance = 1e-10;
    opts.method.form = 'nonlinearprogram';
    opts.method.derivatives = 'complex';

    case 3
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 400; % number of nodes
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.solver.tolerance = 1e-10;
    opts.method.form = 'nonlinearprogram';
    opts.method.derivatives = 'complex';

end

end