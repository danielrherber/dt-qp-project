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
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = LQRScalarTransfer(varargin)

% set p and opts (see LQRScalarTransfer_opts)
% input arguments can be provided in the format 'LQRScalarTransfer(p,opts)'
[p,opts] = DTQP_standardizedinputs(@LQRScalarTransfer_opts,varargin);

%% tunable parameters
% "Hyper-Sensitive"-like behavior
tf = 10000; % final time, requires high-precision solution
p.a = -1; % state matrix
p.b = 1; % input matrix
p.c = 1.5; % initial state
p.d = 1; % final state
p.q = 1; % quadratic state cost
p.r = 1; % quadratic control cost

% % "Energy-Optimal Control"
% tf = 5; % final time, tf > 0
% p.a = 2; % state matrix, a > 0
% p.b = 3; % input matrix, b > 0
% p.c = 1; % initial state, c > 0
% p.d = 5; % final state, d < exp(a*tf)*c
% p.q = 0; % quadratic state cost
% p.r = 1/2; % quadratic control cost

% % "Mass-Spring" Problem
% tf = pi/2; % final time
% p.a = 0; % state matrix
% p.b = 1; % input matrix
% p.c = 0; % initial state
% p.d = 1; % final state
% p.q = -1; % quadratic state cost
% p.r = 1; % quadratic control cost

%% setup
t0 = 0;

% system dynamics
A = p.a;
B = p.b;

% Lagrange term
L(1).left = 1; % control variables
L(1).right = 1; % control variables
L(1).matrix = p.r;
L(2).left = 2; % state variables
L(2).right = 2; % state variables
L(2).matrix = p.q;

% initial conditions
LB(1).right = 4; % initial states
LB(1).matrix = p.c;
UB(1).right = 4; % initial states
UB(1).matrix = p.c;

% final conditions
LB(2).right = 5; % final states
LB(2).matrix = p.d;
UB(2).right = 5; % final states
UB(2).matrix = p.d;

% combine
setup.A = A; setup.B = B; setup.L = L;
setup.LB = LB; setup.UB = UB; setup.t0 = t0; setup.tf = tf; setup.p = p;

%% solve
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = LQRScalarTransfer_output(T,U,Y,P,F,in,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
LQRScalarTransfer_plot(T,U,Y,P,F,in,opts,sol)

end
% User options function for LQRScalarTransfer example
function opts = LQRScalarTransfer_opts
% test number
num = 1;

switch num
case 1
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'CGL';
    opts.dt.nt = 2000; % number of nodes
end
end