%--------------------------------------------------------------------------
% LQRScalar.m
% Scalar linear quadratic regulator
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = LQRScalar(varargin)
% input arguments can be provided in the format 'LQRScalar(auxdata,opts)'

% set local functions
ex_opts = @LQRScalar_opts; % options function
ex_output = @LQRScalar_output; % output function
ex_plot = @LQRScalar_plot; % plot function

% set auxdata and opts (see local_opts)
[auxdata,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
t0 = 0; tf = 1; % time horizon
auxdata.x0 = 10; % initial state
auxdata.a = -1;
auxdata.b = 1;
auxdata.q = 1; % q > -a^2*r/b^2
auxdata.r = 11;
auxdata.m = 10;

%% setup
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

% Mayer term
M(1).left = 5; % final states
M(1).right = 5; % final states
M(1).matrix = auxdata.m;

% initial conditions
LB(1).right = 4; % initial states
LB(1).matrix = auxdata.x0;
UB(1).right = 4; % initial states
UB(1).matrix = auxdata.x0;

% combine
setup.A = A; setup.B = B; setup.L = L; setup.M = M;
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
function opts = LQRScalar_opts
% test number
num = 1;

switch num
case 1
    opts = [];

end

end