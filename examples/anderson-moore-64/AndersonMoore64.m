%--------------------------------------------------------------------------
% AndersonMoore64.m
% p. 64 of B. D. O. Anderson and J. B. Moore, Optimal Control: Linear
% Quadratic Methods. Prentice-Hall, 1989, isbn: 0136386512
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = AndersonMoore64(varargin)

% set p and opts (see AndersonMoore64_opts)
% input arguments can be provided in the format 'AndersonMoore64(p,opts)'
[p,opts] = DTQP_standardizedinputs(@AndersonMoore64_opts,varargin);

%% tunable parameters
t0 = 0; tf = 10; % time horizon
p.x0 = 5; % initial state

%% setup
% system dynamics
A = 1/2; 
B = 1;

% Lagrange term
L(1).left = 1; % control variables
L(1).right = 1; % control variables
L(1).matrix = {@(t) 2*exp(-t)};
L(2).left = 2; % state variables
L(2).right = 2; % state variables
L(2).matrix = {@(t) 0.5*exp(-t)};

% initial states
LB(1).right = 4; % initial states
LB(1).matrix = p.x0;
UB(1).right = 4; % initial states
UB(1).matrix = p.x0;

% combine structures
setup.A = A; setup.B = B; setup.L = L;
setup.LB = LB; setup.UB = UB; setup.t0 = t0; setup.tf = tf; setup.p = p;

%% solve
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = AndersonMoore64_output(T,U,Y,P,F,in,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
AndersonMoore64_plot(T,U,Y,P,F,in,opts,sol)

end
% User options function for AndersonMoore64 example
function opts = AndersonMoore64_opts
% test number
num = 1;

switch num
case 1
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 30;
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