%--------------------------------------------------------------------------
% Tumor.m
% U. Ledzewicz and H. Schättler, Analysis of optimal controls for a
% mathematical model of tumour anti‐angiogenesis, Optim. Control Appl.
% Meth., vol. 29, pp. 41-57, 2008, doi:10.1002/oca.814
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
A0 = 15;
b = 5.85; % per day
Mew = 0.02; % per day
G = 0.15; % per mg of dose per day
zeta = 0.084; % per day
D = 0.00873; % per mm^2 per day
p0 = (((b-Mew)/D)^(3/2))/2;
q0 = p0/2;
umax = 75;
y0 = [p0,q0,0];
ymin = [0.1,0.1,-inf];
auxdata.y0 = y0;
setup.auxdata = auxdata;

% time horizon
setup.t0 = 0; setup.tf = 1.2;

% counts for controls, states, and parameters
setup.counts.nu = 1;
setup.counts.nx = 3;

% Mayer term
M(1).right = 5; % final states
M(1).left = 0; % singleton
M(1).matrix = [1,0,0];
setup.lq.mayer = M;

% system dynamics
setup.nonlin.dynamics = '[-zeta*y1*log(y1/y2); y2*(b-(Mew+(D*(y1^(2/3)))+G*u1)); u1]';

% symbolic data for nonlin
setup.nonlin.data.symbols = 'zeta b Mew D G';
setup.nonlin.data.values = [zeta b Mew D G];

% simple bounds
UB(1).right = 4; UB(1).matrix = y0'; % initial states
LB(1).right = 4; LB(1).matrix = y0';
UB(2).right = 1; UB(2).matrix = umax; % controls
LB(2).right = 1; LB(2).matrix = 0;
UB(3).right = 5; UB(3).matrix = [inf,inf,A0]'; % final states
LB(3).right = 2; LB(3).matrix = ymin'; % states
setup.lq.ub = UB; setup.lq.lb = LB;

% guess
Y0 = [y0;y0];
U0 = [umax;umax];
setup.method.guess.X = [U0,Y0];

% scaling
scaling(1).right = 1; % controls
scaling(1).matrix = umax;
scaling(2).right = 2; % states
scaling(2).matrix = [p0,q0,A0];
setup.method.scaling = scaling;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = Tumor_output(T,U,Y,P,F,in,opts);

% plots
Tumor_plot(T,U,Y,P,F,in,opts,sol)

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
    opts.dt.nt = 1000; % number of nodes
    opts.method.sqpflag = false;
    opts.method.trustregionflag = true;
    opts.method.improveguess = false;
    opts.method.delta = 3000;
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.solver.maxiters = 4000;
case 2
    opts.general.displevel = 1;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100; % number of nodes
    opts.method.form = 'qlin';
    opts.solver.display = 'none';
    opts.method.trustregionflag = false;
    opts.method.improveguess = false;
end

end