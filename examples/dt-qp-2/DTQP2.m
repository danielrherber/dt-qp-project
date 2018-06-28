%--------------------------------------------------------------------------
% DTQP2.m
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = DTQP2(varargin)

% set p and opts (see DTQP2_opts)
% input arguments can be provided in the format 'DTQP2(p,opts)'
[p,opts] = DTQP_standardizedinputs(@DTQP2_opts,varargin);

%% tunable parameters
tf = 15; % time horizon
p.x0 = 1; % initial state
p.r = 1;
p.m = 1;
p.a = 1;
p.b = 1;
p.omega = pi;

%% setup
t0 = 0;

% system dynamics
A = p.a;
B = @(t) p.b*sin(p.omega*t);

% Lagrange term
L(1).left = 1; % control variables
L(1).right = 1; % control variables
L(1).matrix = p.r;

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
setup.LB = LB; setup.UB = UB; setup.t0 = t0; setup.tf = tf; setup.p = p;

%% solve
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = DTQP2_output(T,U,Y,P,F,in,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
DTQP2_plot(T,U,Y,P,F,in,opts,sol)

end
% User options function for DTQP2 example
function opts = DTQP2_opts
% test number
num = 1;

switch num
case 1
    opts.dt.nt = 1000; % number of time points
end
end