%--------------------------------------------------------------------------
% Hypersensitive.m
% A. V. Rao and K. D. Mease, "Eigenvector approximate dichotomic basis
% method for solving hyper-sensitive optimal control problems", *Optimal
% Control Applications and Methods*, vol. 21, no. 1, pp. 1-19, 2000,
% doi: 10.1002/(SICI)1099-1514(200001/02)21:1<1::AID-OCA646>3.0.CO;2-V
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = Hypersensitive(varargin)
% input arguments can be provided in the format 'Hypersensitive(auxdata,opts)'

% set local functions
ex_opts = @Hypersensitive_opts; % options function
ex_output = @Hypersensitive_output; % output function
ex_plot = @Hypersensitive_plot; % plot function

% set auxdata and opts (see local_opts)
[auxdata,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters

%% setup
% time horizon
auxdata.t0 = 0; auxdata.tf = 10000;

% number of controls, states, and parameters
n.nu = 1; n.ny = 1;

% system dynamics
element.dynamics = '[-y1^3+u1]';

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = 1/2;
L(2).left = 2; L(2).right = 2; L(2).matrix = 1/2;

% initial value constraints
UB(1).right = 4; UB(1).matrix = 1.5; % initial states
LB(1).right = 4; LB(1).matrix = 1.5;

% final value constraints
UB(2).right = 5; UB(2).matrix = 1; % final states
LB(2).right = 5; LB(2).matrix = 1;

% combine structures
setup.element = element; setup.L = L; setup.UB = UB; setup.LB = LB;
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
function opts = Hypersensitive_opts
% test number
num = 1;

switch num
case 1 % ipfmincon method
    opts.general.displevel = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    S2 = 10; N2 = 300; p.t0 = 0; p.tf = 10000;
    opts.dt.t = unique([linspace(p.t0,p.t0+S2,N2),linspace(p.t0+S2,p.tf-S2,N2),linspace(p.tf-S2,p.tf,N2)]);
    opts.dt.mesh = 'USER';
    opts.dt.nt = length(opts.dt.t);
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
case 2 % qlin method
    opts.general.displevel = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    S2 = 10; N2 = 300; p.t0 = 0; p.tf = 10000;
    opts.dt.t = unique([linspace(p.t0,p.t0+S2,N2),linspace(p.t0+S2,p.tf-S2,N2),linspace(p.tf-S2,p.tf,N2)]);
    opts.dt.mesh = 'USER';
    opts.dt.nt = length(opts.dt.t);
    opts.solver.display = 'none';
    opts.method.form = 'qlin';
    opts.method.trustregion = false;
    opts.method.sqpflag = false;
end

end