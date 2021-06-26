%--------------------------------------------------------------------------
% DTQP2.m
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = DTQP2(varargin)
% input arguments can be provided in the format 'DTQP2(auxdata,opts)'

% set local functions
ex_opts = @DTQP2_opts; % options function
ex_output = @DTQP2_output; % output function
ex_plot = @DTQP2_plot; % plot function

% set auxdata and opts (see local_opts)
[auxdata,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
tf = 15; % time horizon
auxdata.x0 = 1; % initial state
auxdata.r = 1;
auxdata.m = 1;
auxdata.a = 1;
auxdata.b = 1;
auxdata.omega = pi;

%% setup
t0 = 0;

% system dynamics
A = auxdata.a;
B = @(t) auxdata.b*sin(auxdata.omega*t);

% Lagrange term
L(1).left = 1; % control variables
L(1).right = 1; % control variables
L(1).matrix = auxdata.r;

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
function opts = DTQP2_opts
% test number
num = 2;

switch num
case 1
    opts.dt.nt = 1000; % number of time points
case 2
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 10;
    opts.solver.display = 'none';
    opts.dt.meshr.method = 'SS-BETTS';
    opts.dt.meshr.tolerance = 1e-5;
end

end