%--------------------------------------------------------------------------
% NeuenhofenKerriganX2.m
% X2 from M. P. Neuenhofen and E. C. Kerrigan, "An integral penalty-barrier
% direct transcription method for optimal control," 2020, arXiv:2009.06217
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = NeuenhofenKerriganX2(varargin)
% input arguments can be provided in the format 'NeuenhofenKerriganX2(p,opts)'

% set local functions
ex_opts = @NeuenhofenKerriganX2_opts; % options function
ex_output = @NeuenhofenKerriganX2_output; % output function
ex_plot = @NeuenhofenKerriganX2_plot; % plot function

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters

%% setup
% time horizon
p.t0 = 0; p.tf = 1;

% number of controls, states, and parameters
n.nu = 1; n.ny = 2;

% system dynamics
element.dynamics = '[u1/(2*y1); 4*y1^4 + u1^2]';

% Mayer term
M(1).left = 0; M(1).right = 5; M(1).matrix = [0 1]; % states

% initial value constraints
UB(1).right = 4; UB(1).matrix = [1 0]; % initial states
LB(1).right = 4; LB(1).matrix = [1 0];

% control bounds
LB(2).right = 1; LB(2).matrix = -1; % controls

% state bounds
LB(3).right = 2; LB(3).matrix = [sqrt(0.4),-inf]; % states

% combine structures
setup.element = element; setup.M = M; setup.UB = UB; setup.LB = LB;
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
function opts = NeuenhofenKerriganX2_opts
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
    opts.dt.nt = 200; % number of nodes
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
    opts.dt.nt = 50; % number of nodes
    opts.solver.tolerance = 1e-12;
    opts.solver.maxiters = 200;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
end

end