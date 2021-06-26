%--------------------------------------------------------------------------
% ChemicalReactor.m
% S. J. Citron, Elements of Optimal Control, Holt, Rinehart and Winston,
% New York, 1969
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = ChemicalReactor(varargin)
% input arguments can be provided in the format 'ChemicalReactor(auxdata,opts)'

% set local functions
ex_opts = @ChemicalReactor_opts; % options function
ex_output = @ChemicalReactor_output; % output function
ex_plot = @ChemicalReactor_plot; % plot function

% set auxdata and opts (see local_opts)
[auxdata,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
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

rho = 2.5;
y0 = [1;0.01];

%% setup
% time horizon
auxdata.t0 = 0; auxdata.tf = tf;

% number of controls, states, and parameters
n.nu = 1; n.ny = 2;

% Mayer term
M(1).left = 0; % singleton
M(1).right = 5; % final states
M(1).matrix = [0,-1];

% system dynamics
element.dynamics = '[-u1*y1; u1*y1 - rho*u1^kc*y2]';
element.parameter_list = 'rho kc';
element.parameter_values = [rho kc];

% simple bounds
UB(1).right = 4; UB(1).matrix = y0'; % initial states
LB(1).right = 4; LB(1).matrix = y0';
UB(2).right = 1; UB(2).matrix = au; % controls
LB(2).right = 1; LB(2).matrix = al;
UB(3).right = 2; UB(3).matrix = [1.1;1.1]; % states
LB(3).right = 2; LB(3).matrix = [-0.1;-0.1];

% guess
Y0 = [[y0'];[y0']];
U0 = [[au];[al]];
setup.guess.X = [U0,Y0];

% scaling
setup.scaling(1).right = 1; % controls
setup.scaling(1).matrix = [au];
setup.scaling(2).right = 2; % states
setup.scaling(2).matrix = [1.1,1.1];

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
function opts = ChemicalReactor_opts
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