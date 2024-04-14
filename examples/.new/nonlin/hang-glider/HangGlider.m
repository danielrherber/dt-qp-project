%--------------------------------------------------------------------------
% HangGlider.m
% R. Bulirsch, E. Nerz, H. J. Pesch, and O. von Stryk, "Combining Direct
% and Indirect Methods in Optimal Control: Range Maximization of a Hang
% Glider,*" in Optimal Control, Birkhäuser Basel, 1993, pp. 273–288
% [Online]. Available: http://dx.doi.org/10.1007/978-3-0348-7539-4_20
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
rho = 1.13;
CD0 = 0.034;
k = 0.069662;
g = 9.80665;
m = 100;
S = 14;
uM = 2.5;
R = 100;
umax = 1.4;
y0 = [0,1000,13.227567500,-1.2875005200];
yf = [inf,900,13.227567500,-1.2875005200];

% time horizon
setup.t0 = 0; setup.tf = 1;

% number of controls, states, and parameters
setup.counts.nu = 1;
setup.counts.nx = 4;
setup.counts.np = 1;

% Mayer term
M(1).right = 5; % final states
M(1).left = 0; % singleton
M(1).matrix = [-1,0,0,0];
setup.lq.mayer = M;

% system dynamics
strD{1} = '[';
strD{end+1} = 'p1*y3;';
strD{end+1} = 'p1*y4;';
strD{end+1} = 'p1*( - (S*rho*y3*(k*u1^2 + CD0)*((y4 + uM*exp(-(y1/R - 5/2)^2)*((y1/R - 5/2)^2 - 1))^2 + y3^2)^(1/2))/(2*m) - (S*rho*u1*(y4 + uM*exp(-(y1/R - 5/2)^2)*((y1/R - 5/2)^2 - 1))*((y4 + uM*exp(-(y1/R - 5/2)^2)*((y1/R - 5/2)^2 - 1))^2 + y3^2)^(1/2))/(2*m) );';
strD{end+1} = 'p1*( (S*rho*u1*y3*((y4 + uM*exp(-(y1/R - 5/2)^2)*((y1/R - 5/2)^2 - 1))^2 + y3^2)^(1/2))/(2*m) - g - (S*rho*(k*u1^2 + CD0)*(y4 + uM*exp(-(y1/R - 5/2)^2)*((y1/R - 5/2)^2 - 1))*((y4 + uM*exp(-(y1/R - 5/2)^2)*((y1/R - 5/2)^2 - 1))^2 + y3^2)^(1/2))/(2*m) )';
strD{end+1} = ']';
setup.nonlin.dynamics = horzcat(strD{:});

% symbolic data for nonlin
setup.nonlin.data.symbols = 'rho CD0 k g m S uM R';
setup.nonlin.data.values = [rho CD0 k g m S uM R];

% simple bounds
UB(1).right = 4; UB(1).matrix = y0'; % initial states
LB(1).right = 4; LB(1).matrix = y0';
UB(2).right = 1; UB(2).matrix = umax; % controls
LB(2).right = 1; LB(2).matrix = 0;
UB(3).right = 5; UB(3).matrix = [inf,yf(2:4)]'; % final states
LB(3).right = 5; LB(3).matrix = [-inf,yf(2:4)]';
UB(4).right = 2; UB(4).matrix = [3000 3000 15 15]; % states
LB(4).right = 2; LB(4).matrix = [-3000 0 -15 -15];
UB(5).right = 3; UB(5).matrix = 200; % parameters
LB(5).right = 3; LB(5).matrix = 0;
setup.lq.ub = UB; setup.lq.lb = LB;

% guess
Y0 = [y0;[1250,yf(2:4)]];
U0 = [1;1];
P0 = [100;100];
setup.method.guess.X = [U0,Y0,P0];

% scaling
scaling(1).right = 1; % controls
scaling(1).matrix = umax;
scaling(2).right = 2; % states
scaling(2).matrix = [3000 3000 15 15];
scaling(3).right = 3; % parameters
scaling(3).matrix = 200;
setup.method.scaling = scaling;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = HangGlider_output(T,U,Y,P,F,in,opts);

% plots
HangGlider_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts

% test number
num = 4;

switch num
case 1
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 500; % number of nodes
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.solver.maxiters = 20000;
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 200; % number of nodes
    opts.method.sqpflag = false;
    opts.method.trustregionflag = true;
    opts.method.improveguess = false;
    opts.method.delta = 3000;
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.solver.maxiters = 20000;
case 3
    opts.general.displevel = 1;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 1000; % number of nodes
    opts.method.form = 'qlin';
    opts.solver.display = 'none';
    opts.method.trustregionflag = false;
    opts.method.improveguess = false;
case 4
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100; % number of nodes
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.solver.maxiters = 20000;
    opts.solver.tolerance = 1e-12;
    opts.dt.meshr.method = 'RICHARDSON-DOUBLING';
    opts.dt.meshr.tolerance = 1e-6;
end

end