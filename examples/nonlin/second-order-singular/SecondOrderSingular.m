%--------------------------------------------------------------------------
% SecondOrderSingular.m
% G. M. Aly, "The computation of optimal singular control", International
% Journal of Control, vol. 28, no. 5, pp. 681–688, Nov. 1978,
% doi: 10.1080/00207177808922489
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = SecondOrderSingular(varargin)
% input arguments can be provided in the format 'SecondOrderSingular(auxdata,opts)'

% set local functions
ex_opts = @SecondOrderSingular_opts;
ex_output = @SecondOrderSingular_output;
ex_plot = @SecondOrderSingular_plot;

% set auxdata and opts (see local_opts)
[auxdata,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
tf = 5;

%% setup
% time horizon
auxdata.t0 = 0; auxdata.tf = tf;

% number of controls, states, and parameters
n.nu = 1; n.ny = 3;

% system dynamics
str{1} = '[';
str{end+1} = 'y2; ';
str{end+1} = 'u1; ';
str{end+1} = 'y1^2/2 + y2^2/2';
str{end+1} = ']';
element.dynamics = horzcat(str{:});

% Mayer term
M(1).left = 0; M(1).right = 5; M(1).matrix = [0 0 1]; % final states

% simple bounds
UB(1).right = 4; UB(1).matrix = [0,1,0]; % initial states
LB(1).right = 4; LB(1).matrix = [0,1,0];
UB(2).right = 1; UB(2).matrix = 1; % controls
LB(2).right = 1; LB(2).matrix = -1;

% guess
Y0 = [[0,1,0];[0,1,1]];
U0 = [[-1];[-1]];
setup.guess.X = [U0,Y0];

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
function opts = SecondOrderSingular_opts
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
    opts.solver.maxiters = 300;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.method.olqflag = true;
    opts.method.derivatives = 'complex';
    opts.solver.tolerance = 1e-12;
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 100; % number of nodes
    opts.solver.maxiters = 300;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.method.olqflag = true;
    opts.method.derivatives = 'complex';
end

end