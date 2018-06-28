%--------------------------------------------------------------------------
% DTQP3.m
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = DTQP3(varargin)

% set p and opts (see DTQP3_opts)
% input arguments can be provided in the format 'DTQP3(p,opts)'
[p,opts] = DTQP_standardizedinputs(@DTQP3_opts,varargin);

%% tunable parameters
tf = 10; % time horizon
p.x0 = 1; % initial state
p.b = 1;
p.r = 1;
p.m = 1;
p.a1 = 2;
p.a2 = 1;
p.w1 = 3;
p.w2 = 8;
p.ParameterFlag = 1; % parameter version?

%% setup
t0 = 0;

% system dynamics
if p.ParameterFlag
	setup.G = p.b;
else
	setup.B = p.b;
end
setup.d{1,1} = @(t) p.a1*sin(p.w1*t) + p.a2*sin(p.w2*t);

% Lagrange term
if p.ParameterFlag
    L(1).left = 3; % parameters 
    L(1).right = 3; % parameters
    L(1).matrix = p.r;
else
    L(1).left = 1; % control variables
    L(1).right = 1; % control variables
    L(1).matrix = p.r;
end

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
setup.L = L; setup.M = M; setup.LB = LB; setup.UB = UB;
setup.t0 = t0; setup.tf = tf; setup.p = p;

%% solve
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = DTQP3_output(T,U,Y,P,F,in,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
DTQP3_plot(T,U,Y,P,F,in,opts,sol)

end
% User options function for DTQP3 example
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