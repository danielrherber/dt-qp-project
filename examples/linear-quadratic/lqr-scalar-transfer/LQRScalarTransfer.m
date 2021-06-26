%--------------------------------------------------------------------------
% LQRScalarTransfer.m
%
% Similar behavior to the "Hyper-Sensitive" problem:
% A. V. Rao and K. D. Mease, Eigenvector Approximate Dichotomic Basis
% Methods for Solving Hyper-Sensitive Optimal Control Problems, Optimal
% Control Applications and Methods, Vol. 21, No. 1., January-February,
% 2000, pp. 1-17.
%
% "Energy-Optimal Control" problem in reference below is a special case:
% H. P. Geering, Optimal Control with Engineering Applications, Springer,
% 2007, pp. 46-48, doi: 10.1007/978-3-540-69438-0
%
% "Mass-Spring" problem included in GPOPS-II is a special case.
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = LQRScalarTransfer(varargin)
% input arguments can be provided in the format 'LQRScalarTransfer(auxdata,opts)'

% set local functions
ex_opts = @LQRScalarTransfer_opts; % options function
ex_output = @LQRScalarTransfer_output; % output function
ex_plot = @LQRScalarTransfer_plot; % plot function

% set auxdata and opts (see local_opts)
[auxdata,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
% "Hyper-Sensitive"-like behavior
tf = 10000; % final time, requires high-precision solution
auxdata.a = -1; % state matrix
auxdata.b = 1; % input matrix
auxdata.c = 1.5; % initial state
auxdata.d = 1; % final state
auxdata.q = 1; % quadratic state cost
auxdata.r = 1; % quadratic control cost

% % "Energy-Optimal Control"
% tf = 5; % final time, tf > 0
% auxdata.a = 2; % state matrix, a > 0
% auxdata.b = 3; % input matrix, b > 0
% auxdata.c = 1; % initial state, c > 0
% auxdata.d = 5; % final state, d < exp(a*tf)*c
% auxdata.q = 0; % quadratic state cost
% auxdata.r = 1/2; % quadratic control cost

% % "Mass-Spring" Problem
% tf = pi/2; % final time
% auxdata.a = 0; % state matrix
% auxdata.b = 1; % input matrix
% auxdata.c = 0; % initial state
% auxdata.d = 1; % final state
% auxdata.q = -1; % quadratic state cost
% auxdata.r = 1; % quadratic control cost

%% setup
t0 = 0;

% system dynamics
A = auxdata.a;
B = auxdata.b;

% Lagrange term
L(1).left = 1; % control variables
L(1).right = 1; % control variables
L(1).matrix = auxdata.r;
L(2).left = 2; % state variables
L(2).right = 2; % state variables
L(2).matrix = auxdata.q;

% initial conditions
LB(1).right = 4; % initial states
LB(1).matrix = auxdata.c;
UB(1).right = 4; % initial states
UB(1).matrix = auxdata.c;

% final conditions
LB(2).right = 5; % final states
LB(2).matrix = auxdata.d;
UB(2).right = 5; % final states
UB(2).matrix = auxdata.d;

% combine
setup.A = A; setup.B = B; setup.L = L;
setup.LB = LB; setup.UB = UB; setup.t0 = t0; setup.tf = tf; setup.auxdata = auxdata;

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
function opts = LQRScalarTransfer_opts
% test number
num = 2;

switch num
case 1
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'CGL';
    opts.dt.nt = 2000; % number of nodes
case 2
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 10; % number of nodes
    opts.solver.tolerance = 1e-15;
    opts.solver.display = 'none';
    opts.dt.meshr.method = 'SS-BETTS';
    opts.dt.meshr.tolerance = 1e-6;
    opts.dt.meshr.ntmaxinterval = 5;
end

end