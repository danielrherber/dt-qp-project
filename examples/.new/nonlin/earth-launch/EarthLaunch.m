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

% auxiliary data
auxdata.t0 = 0; auxdata.tf = 1;
setup.auxdata = auxdata;

% time horizon
setup.t0 = auxdata.t0; setup.tf = auxdata.tf;

% number of controls, states, and parameters
setup.counts.nu = 1;
setup.counts.nx = 5;
setup.counts.np = 1;

% Mayer term
M(1).right = 3; % parameters
M(1).left = 0; % singleton
M(1).matrix = 1;
setup.lq.mayer = M;

% system dynamics
strD{1} = '[';
strD{end+1} = 'p1*(y3);'; % x, range [m]
strD{end+1} = 'p1*(y4);'; % y, altitude [m]
strD{end+1} = 'p1*((FT*cos(u1) - 0.5*CDA*(rho_ref*exp(-y2/hs))*y3^2)/y5);'; % vx, x velocity [m/s]
strD{end+1} = 'p1*((FT*sin(u1) - 0.5*CDA*(rho_ref*exp(-y2/hs))*y4^2)/y5 - g);'; % vy, y velocity [m/s]
strD{end+1} = 'p1*(-FT/(g*Isp));'; % m, mass [kg]
strD{end+1} = ']';
setup.nonlin.dynamics = horzcat(strD{:});

% data for symbolic functions
setup.nonlin.data.symbols = 'FT CDA g Isp rho_ref hs';
setup.nonlin.data.values = [FT CD*S g Isp rho_ref h_ref];

% simple bounds
UB(1).right = 4; UB(1).matrix = [x0,y0,vx0,vy0,m0]'; % initial states
LB(1).right = 4; LB(1).matrix = [x0,y0,vx0,vy0,m0]';
UB(2).right = 1; UB(2).matrix = pi/2; % controls
LB(2).right = 1; LB(2).matrix = -pi/2;
UB(3).right = 5; UB(3).matrix = [inf,yf,vxf,vyf,inf]'; % final states
LB(3).right = 5; LB(3).matrix = [-inf,yf,vxf,vyf,-inf]';
UB(4).right = 3; UB(4).matrix = 500; % parameters
LB(4).right = 3; LB(4).matrix = 10;
setup.lq.ub = UB; setup.lq.lb = LB;

% guess
Y0 = [[x0,y0,vx0,vy0,m0];[1.15e5,yf,vxf,vyf,1163]];
U0 = [1.5;-0.76];
P0 = [150;150];
setup.method.guess.X = [U0,Y0,P0];

% scaling
scaling(1).right = 1; % controls
scaling(1).matrix = pi/2;
scaling(2).right = 2; % states
scaling(2).matrix = [hs,hs,vs,vs,m0];
scaling(3).right = 1; % parameters
scaling(3).matrix = ts;
setup.method.scaling = scaling;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = EarthLaunch_output(T,U,Y,P,F,in,opts);

% plots
EarthLaunch_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts

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