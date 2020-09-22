%--------------------------------------------------------------------------
% NeuenhofenKerriganX1.m
% X1 from M. P. Neuenhofen and E. C. Kerrigan, "An integral penalty-barrier
% direct transcription method for optimal control," 2020, arXiv:2009.06217
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = NeuenhofenKerriganX1(varargin)
% input arguments can be provided in the format 'NeuenhofenKerriganX1(p,opts)'

% set local functions
ex_opts = @NeuenhofenKerriganX1_opts; % options function
ex_output = @NeuenhofenKerriganX1_output; % output function
ex_plot = @NeuenhofenKerriganX1_plot; % plot function

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters

%% setup
% time horizon
p.t0 = 0; p.tf = pi/2;

% number of controls, states, and parameters
n.nu = 1; n.ny = 1;

% system dynamics
symb.D = '[y1^2/2+u1]';

% Lagrange term
L(1).left = 0; L(1).right = 1; L(1).matrix = @(t) cos(t); % controls
L(2).left = 2; L(2).right = 2; L(2).matrix = 1; % states

% initial value constraints
UB(1).right = 4; UB(1).matrix = 0; % initial states
LB(1).right = 4; LB(1).matrix = 0;

% control bounds
UB(2).right = 1; UB(2).matrix = 1; % controls
LB(2).right = 1; LB(2).matrix = -1.5;

% combine structures
setup.symb = symb; setup.L = L; setup.UB = UB; setup.LB = LB;
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
function opts = NeuenhofenKerriganX1_opts
% test number
num = 2;

switch num
case 1
    opts.general.plotflag = 1; % create the plots
    opts.general.saveflag = false;
    opts.general.displevel = 2;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100; % number of nodes
    opts.solver.tolerance = 1e-12;
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
    opts.solver.tolerance = 1e-12;
    opts.solver.maxiters = 200;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
end

end