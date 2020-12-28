%--------------------------------------------------------------------------
% Nonlinear1D.m
% D. Garg et al., "Direct Trajectory Optimization and Costate Estimation of
% General Optimal Control Problems Using a Radau Pseudospectral Method," in
% AIAA Guidance, Navigation, and Control Conference, 2009,
% doi: 10.2514/6.2009-5989
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = Nonlinear1D(varargin)
% input arguments can be provided in the format 'Nonlinear1D(p,opts)'

% set local functions
ex_opts = @Nonlinear1D_opts; % options function
ex_output = @Nonlinear1D_output; % output function
ex_plot = @Nonlinear1D_plot; % plot function

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters

%% setup
% time horizon
p.t0 = 0; p.tf = 5;

% number of controls, states, and parameters
n.nu = 1; n.ny = 1;

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = 1/2;
L(2).left = 0; L(2).right = 2; L(2).matrix = 1/2;

% system dynamics
symb.D = '[2*y1 + 2*u1*sqrt(y1)]';

% simple bounds
UB(1).right = 4; UB(1).matrix = [2]; % initial states
LB(1).right = 4; LB(1).matrix = [2];
UB(2).right = 5; UB(2).matrix = [1]; % final states
LB(2).right = 5; LB(2).matrix = [1];
LB(3).right = 2; LB(3).matrix = [0]; % optional constraint for sqrt(y1)

% guess
Y0 = [[2];[1]];
U0 = [[0];[0]];
p.guess = [U0,Y0];

% combine structures
setup.symb = symb; setup.UB = UB; setup.LB = LB; setup.L = L;
setup.t0 = p.t0; setup.tf = p.tf; setup.p = p; setup.n = n;

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
function opts = Nonlinear1D_opts
% test number
num = 2;

switch num
case 1
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100; % number of nodes
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.solver.tolerance = 1e-10;
    opts.dt.meshr.method = 'RICHARDSON-DOUBLING';
    opts.dt.meshr.tolerance = 1e-6;
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 8; % number of nodes
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.solver.tolerance = 1e-10;
    opts.dt.meshr.method = 'RICHARDSON-DOUBLING';
    opts.dt.meshr.tolerance = 1e-6;
end

end