%--------------------------------------------------------------------------
% SimpleSASA.m
% D. R. Herber and J. T. Allison, "Unified Scaling of Dynamic Optimization
% Design Formulations," in Volume 2A: 43rd Design Automation Conference,
% 2017, doi: 10.1115/detc2017-67676
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = SimpleSASA(varargin)
% input arguments can be provided in the format 'SimpleSASA(p,opts)'

% set local functions
ex_opts = @SimpleSASA_opts; % options function
ex_output = @SimpleSASA_output; % output function
ex_plot = @SimpleSASA_plot; % plot function

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
p.umax = 1;
p.J = 1;
tf = 2;

%% setup
% time horizon
p.t0 = 0; p.tf = tf;

% number of controls, states, and parameters
n.nu = 1; n.ny = 2; n.np = 1;

% system dynamics
symb.D = '[y2; -p1/J*y1 + u1/J]';
symb.paramstr = 'J';
symb.param = [p.J];

% Mayer term
M(1).right = 5; M(1).left = 0; M(1).matrix = [-1,0]; % final states

% simple bounds
UB(1).right = 4; UB(1).matrix = [0,0]; % initial states
LB(1).right = 4; LB(1).matrix = [0,0];
UB(2).right = 5; UB(2).matrix = [inf 0]; % final states
LB(2).right = 5; LB(2).matrix = [-inf 0];
UB(3).right = 1; UB(3).matrix = p.umax; % controls
LB(3).right = 1; LB(3).matrix = -p.umax;

% guess
Y0 = [[0,0];[1 0]];
U0 = [[0];[0]];
P0 = [[1];[1]];
p.guess = [U0,Y0,P0];

% combine structures
setup.symb = symb; setup.M = M; setup.UB = UB; setup.LB = LB;
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
function opts = SimpleSASA_opts
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
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
end

end