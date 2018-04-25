%--------------------------------------------------------------------------
% LQRScalar.m
% Scalar linear quadratic regular
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = LQRScalar(varargin)

% set p and opts (see LQRScalar_opts.m)
% input arguments can be provided in the format 'LQRScalar(p,opts)'
[p,opts] = DTQP_standardizedinputs('LQRScalar_opts',varargin);

%% tunable parameters
p.t0 = 0; p.tf = 1; % time horizon
p.x0 = 10; % initial state
p.a = -1;
p.b = 1;
p.q = 1; % q > -a^2*r/b^2
p.r = 11;
p.m = 10;

%% setup
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

% Mayer term
M(1).left = 5; % final states
M(1).right = 5; % final states
M(1).matrix = p.m;

% initial conditions
LB(1).right = 4; % initial states
LB(1).matrix = p.x0;
UB(1).right = 4; % initial states
UB(1).matrix = p.x0;

% combine
setup.A = A; setup.B = B; setup.L = L; setup.M = M; 
setup.LB = LB; setup.UB = UB; setup.p = p;

%% solve
[T,U,Y,P,F,p,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = LQRScalar_output(T,U,Y,P,F,p,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
LQRScalar_plot(T,U,Y,P,F,p,opts,sol)