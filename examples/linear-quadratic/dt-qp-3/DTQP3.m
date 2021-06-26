%--------------------------------------------------------------------------
% DTQP3.m
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = DTQP3(varargin)
% input arguments can be provided in the format 'DTQP3(auxdata,opts)'

% set local functions
ex_opts = @DTQP3_opts; % options function
ex_output = @DTQP3_output; % output function
ex_plot = @DTQP3_plot; % plot function

% set auxdata and opts (see local_opts)
[auxdata,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
tf = 10; % time horizon
auxdata.x0 = 1; % initial state
auxdata.b = 1;
auxdata.r = 1;
auxdata.m = 1;
auxdata.a1 = 2;
auxdata.a2 = 1;
auxdata.w1 = 3;
auxdata.w2 = 8;
auxdata.ParameterFlag = 1; % parameter version?

%% setup
t0 = 0;

% system dynamics
if auxdata.ParameterFlag
	setup.G = auxdata.b;
else
	setup.B = auxdata.b;
end
setup.d{1,1} = @(t) auxdata.a1*sin(auxdata.w1*t) + auxdata.a2*sin(auxdata.w2*t);

% Lagrange term
if auxdata.ParameterFlag
    L(1).left = 3; % parameters
    L(1).right = 3; % parameters
    L(1).matrix = auxdata.r;
else
    L(1).left = 1; % control variables
    L(1).right = 1; % control variables
    L(1).matrix = auxdata.r;
end

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
setup.L = L; setup.M = M; setup.LB = LB; setup.UB = UB;
setup.t0 = t0; setup.tf = tf; setup.auxdata = auxdata;

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
function opts = DTQP3_opts
% test number
num = 1;

switch num
case 1
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 100;
case 2
    opts.dt.quadrature = 'CQHS';
    opts.dt.defects = 'HS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100;
case 3
    opts.dt.quadrature = 'CTR';
    opts.dt.defects = 'TR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100;
end

end