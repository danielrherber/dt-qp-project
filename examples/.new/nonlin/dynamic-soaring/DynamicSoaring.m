%--------------------------------------------------------------------------
% DynamicSoaring.m
% Y. J. Zhao, "Optimal Pattern of Glider Dynamic Soaring," Optimal Control
% Applications and Methods, vol. 25, no. 2, pp. 67-89, 2004,
% doi: 10.1002/oca.710
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
rho = 0.002378;
CD0 = 0.00873;
K = 0.045;
g = 32.2;
m = 5.6;
S = 45.09703;
W0 = 0;
lmin = -2;
lmax = 5;

auxdata.rho = rho;
auxdata.S = S;
auxdata.m = m;
auxdata.g = g;
auxdata.lmin = lmin;
auxdata.lmax = lmax;
setup.auxdata = auxdata;

% initial values
x0 = 0;
y0 = 0;
z0 = 0;
v0 = 100;

% variable bounds and scaling factor
tfmin = 1; tfmax = 30; tfscale = tfmax;
xmin = -1000; xmax = 1000; xscale = xmax;
ymin = -1000; ymax = 1000; yscale = ymax;
zmin = 0; zmax = 1000; zscale = zmax;
vmin = 10; vmax = 350;  vscale = vmax;
gammamin = -75*pi/180; gammamax = 75*pi/180; gammascale = gammamax;
psimin = -3*pi; psimax = pi/2; psiscale = abs(psimin);
betamin = 0.005; betamax = 0.15; betascale = betamax;
CLmin = -0.5; CLmax = 1.5; CLscale = CLmax;
Phimin = -75/180*pi; Phimax = 75/180*pi; Phiscale = Phimax;

% time horizon
setup.t0 = 0; setup.tf = 1;

% number of controls, states, and parameters
setup.counts.nu = 2;
setup.counts.nx = 6;
setup.counts.np = 2;
% y1 = x
% y2 = y
% y3 = z
% y4 = v
% y5 = gamma
% y6 = psi
% u1 = CL
% u2 = phi
% p1 = beta
% p2 = tf

% Mayer term
M(1).left = 0; % singleton
M(1).right = 3; % parameters
M(1).matrix = [1 0];
setup.lq.mayer = M;

% system dynamics
str{1} = '[';
str{end+1} = 'p2*(y4*cos(y5)*sin(y6) + (p1*y3 + W0)); '; % x
str{end+1} = 'p2*(y4*cos(y5)*cos(y6)); '; % y
str{end+1} = 'p2*(y4*sin(y5)); '; % z
str{end+1} = 'p2*(-rho*S/2/m*(CD0 + K*u1^2)*y4^2 - (g*1)*sin(y5) - 1*p1*y4*sin(y5)*sin(y6)*cos(y5)); '; % v
str{end+1} = 'p2*(rho*S/2/m*u1*y4*cos(u2) - (g*1)*cos(y5)/y4 + 1*p1*y4*sin(y5)*sin(y6)*sin(y5)/y4); '; % gamma
str{end+1} = 'p2*((rho*S/2/m*u1*y4*sin(u2) - 1*p1*y4*sin(y5)*cos(y6)/y4)/cos(y5))'; % psi
str{end+1} = ']';
setup.nonlin.dynamics = horzcat(str{:});

% data for symbolic functions
setup.nonlin.data.symbols = 'W0 rho S m CD0 K g lmin lmax';
setup.nonlin.data.values = [W0 rho S m CD0 K g lmin lmax];

% nonlinear inequality constraints
setup.nonlin.inequality.func = '[0.5*rho*S/m/g*u1*y4^2 - lmax, lmin - 0.5*rho*S/m/g*u1*y4^2]';
setup.nonlin.inequality.pathboundary = [1 1];

% simple bounds
UB(1).right = 4; UB(1).matrix = [x0, y0, z0, vmax, gammamax, psimax]; % initial states
LB(1).right = 4; LB(1).matrix = [x0, y0, z0, vmin, gammamin, psimin];
UB(2).right = 5; UB(2).matrix = [x0, y0, z0, vmax, gammamax, psimax]; % final states
LB(2).right = 5; LB(2).matrix = [x0, y0, z0, vmin, gammamin, psimin];
UB(3).right = 2; UB(3).matrix = [xmax, ymax, zmax, vmax, gammamax, psimax]; % states
LB(3).right = 2; LB(3).matrix = [xmin, ymin, zmin, vmin, gammamin, psimin];
UB(4).right = 1; UB(4).matrix = [CLmax, Phimax]; % controls
LB(4).right = 1; LB(4).matrix = [CLmin, Phimin];
UB(5).right = 3; UB(5).matrix = [betamax,tfmax]; % parameters
LB(5).right = 3; LB(5).matrix = [betamin,tfmin];
setup.lq.ub = UB; setup.lq.lb = LB;

% linear equality constraints
Y(1).linear(1).right = 5; Y(1).linear(1).matrix = [0;0;0;1;0;0]; % y4(tf)
Y(1).linear(2).right = 4; Y(1).linear(2).matrix = [0;0;0;-1;0;0]; % -y4(t0)
Y(1).b = 0;
Y(2).linear(1).right = 5; Y(2).linear(1).matrix = [0;0;0;0;1;0]; % y5(tf)
Y(2).linear(2).right = 4; Y(2).linear(2).matrix = [0;0;0;0;-1;0]; % -y5(t0)
Y(2).b = 0;
Y(3).linear(1).right = 5; Y(3).linear(1).matrix = [0;0;0;0;0;1]; % y6(tf)
Y(3).linear(2).right = 4; Y(3).linear(2).matrix = [0;0;0;0;0;-1]; % -y6(t0)
Y(3).b = -2*pi;
setup.lq.equality = Y;

% scaling
scaling(1).right = 1; % controls
scaling(1).matrix = [CLscale; Phiscale];
scaling(2).right = 2; % states
scaling(2).matrix = [xscale yscale zscale vscale gammascale psiscale];
scaling(3).right = 3; % parameters
scaling(3).matrix = [betascale; tfscale];
setup.method.scaling = scaling;

% guess
tfguess = 24;
N = 100;
CL0 = CLmax;
tGuess = linspace(0,1,N).';
xguess = 500*cos(2*pi*tGuess)-500;
yguess = 300*sin(2*pi*tGuess);
zguess = -400*cos(2*pi*tGuess)+400;
vguess = 0.8*v0*(1.5+cos(2*pi*tGuess));
gammaguess = pi/6*sin(2*pi*tGuess);
psiguess = -1-tfguess*tGuess/4;
CLguess = CL0*ones(N,1)/3;
phiguess = -ones(N,1);
betaguess = 0.08;

Y0 = [xguess, yguess, zguess, vguess, gammaguess, psiguess];
U0 = [CLguess, phiguess];
P0 = [betaguess*ones(N,1), tfguess*ones(N,1)];
setup.method.guess.X = [U0,Y0,P0];
setup.method.guess.T = tGuess;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = DynamicSoaring_output(T,U,Y,P,F,in,opts);

% plots
DynamicSoaring_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts

% test number
num = 1;

switch num
case 1
    opts.general.plotflag = 1; % create the plots
    opts.general.saveflag = false;
    opts.general.displevel = 2;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 600; % number of nodes
    opts.solver.tolerance = 1e-12;
    opts.solver.maxiters = 500;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.method.derivatives = 'symbolic';
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 50; % number of nodes
    opts.solver.tolerance = 1e-12;
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.solver.display = 'iter';
    opts.solver.maxiters = 2000;
end

end