%--------------------------------------------------------------------------
% ChemicalReactor.m
% S. J. Citron, Elements of Optimal Control, Holt, Rinehart and Winston,
% New York, 1969
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

% different cases
testnum = 1;

switch testnum
    case 1
	al = 0.1; au = 0.5; tf = 2; kc = 1.5;
    case 2
	al = 0.1; au = 0.5; tf = 4; kc = 1.5;
    case 3
	al = 0.1; au = 0.5; tf = 8; kc = 1.5;
    case 4
	al = 0.1; au = 0.2; tf = 2; kc = 1.5;
    case 5
	al = 0.1; au = 0.3; tf = 2; kc = 1.5;
    case 6
	al = 0.1; au = 0.4; tf = 2; kc = 1.5;
    case 7
	al = 0.01; au = 8; tf = 2; kc = 1.5;
    case 8
	al = 0.01; au = 8; tf = 4; kc = 1.5;
    case 9
	al = 0.01; au = 8; tf = 8; kc = 1.5;
    case 10
	al = 0.1; au = 0.5; tf = 2; kc = 0.5;
    otherwise
end

% auxiliary data
rho = 2.5;
y0 = [1;0.01];

% time horizon
setup.t0 = 0; setup.tf = tf;

% number of controls, states, and parameters
setup.counts.nu = 1;
setup.counts.nx = 2;

% Mayer term
M(1).left = 0; % singleton
M(1).right = 5; % final states
M(1).matrix = [0,-1];
setup.lq.mayer = M;

% system dynamics
setup.nonlin.dynamics = '[-u1*y1; u1*y1 - rho*u1^kc*y2]';

% data for symbolic functions
setup.nonlin.data.symbols = 'rho kc';
setup.nonlin.data.values = [rho kc];

% simple bounds
UB(1).right = 4; UB(1).matrix = y0'; % initial states
LB(1).right = 4; LB(1).matrix = y0';
UB(2).right = 1; UB(2).matrix = au; % controls
LB(2).right = 1; LB(2).matrix = al;
UB(3).right = 2; UB(3).matrix = [1.1;1.1]; % states
LB(3).right = 2; LB(3).matrix = [-0.1;-0.1];
setup.lq.ub = UB; setup.lq.lb = LB;

% guess
Y0 = [y0';y0'];
U0 = [au;al];
setup.method.guess.X = [U0,Y0];

% scaling
scaling(1).right = 1; % controls
scaling(1).matrix = au;
scaling(2).right = 2; % states
scaling(2).matrix = [1.1,1.1];
setup.method.scaling = scaling;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = ChemicalReactor_output(T,U,Y,P,F,in,opts);

% plots
ChemicalReactor_plot(T,U,Y,P,F,in,opts,sol)

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
    opts.dt.nt = 100; % number of nodes
    opts.solver.tolerance = 1e-8;
    opts.solver.maxiters = 200;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
case 2
    opts.general.plotflag = 1; % create the plots
    opts.general.saveflag = false;
    opts.general.displevel = 2;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 30; % number of nodes
    opts.solver.tolerance = 1e-8;
    opts.solver.maxiters = 200;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
case 3 % qlin method
    opts.general.plotflag = 1; % create the plots
    opts.general.displevel = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100; % number of nodes
    opts.solver.maxiters = 200;
    opts.solver.display = 'none';
    opts.method.form = 'qlin';
    opts.method.trustregionflag = false;
    opts.method.sqpflag = false;
    opts.method.delta = inf;
    opts.method.improveguess = false; % disabled
end

end