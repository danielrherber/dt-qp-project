%--------------------------------------------------------------------------
% LinearPendulum.m
% Based on the linearized pendulum example in:
% A. Bressan, "Viscosity Solutions of Hamilton-Jacobi Equations and Optimal
% Control Problems," Lecture Notes, pp. 30-31, 2011.
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = LinearPendulum(varargin)
% input arguments can be provided in the format 'LinearPendulum(p,opts)'

% set local functions
ex_opts = @LinearPendulum_opts; % options function
ex_output = @LinearPendulum_output; % output function
ex_plot = @LinearPendulum_plot; % plot function

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
tf = 5;
p.m = 1/6;
p.k = 1;
p.x0 = -4;
p.v0 = 2;
p.umax = 2;

% original parameters
% tf = 20;
% p.m = 1;
% p.k = 1;
% p.x0 = 0;
% p.v0 = 0;
% p.umax = 1;

%% setup
% system dynamics
A = [0 1;-p.k/p.m 0];
B = [0;1/p.m];

% Mayer term
M(1).right = 5; % final states
M(1).matrix = -[1 0];

% initial conditions
LB(1).right = 4; % initial states
LB(1).matrix = [p.x0 p.v0];
UB(1).right = 4; % initial states
UB(1).matrix = [p.x0 p.v0];
LB(2).right = 1; % controls
LB(2).matrix = -p.umax;
UB(2).right = 1; % controls
UB(2).matrix = p.umax;

% combine
setup.A = A; setup.B = B; setup.M = M;
setup.LB = LB; setup.UB = UB; setup.t0 = 0; setup.tf = tf; setup.p = p;

%% solve
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = ex_output(T,U,Y,P,F,in,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
ex_plot(T,U,Y,P,F,in,opts,sol)

end
% User options function for this example
function opts = LinearPendulum_opts
% test number
num = 3;

switch num
case 1
    opts = [];
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 600; % number of nodes
case 2
    opts = [];
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 100; % number of nodes
case 3
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 50;
    opts.solver.display = 'none';
    opts.solver.tolerance = 1e-15;
    opts.dt.meshr.method = 'SS-BETTS';
    opts.dt.meshr.tolerance = 1e-6;
    opts.dt.meshr.ntmaxinterval = 20;
end

end