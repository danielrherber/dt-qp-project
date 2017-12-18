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

% number of time points
p.nt = 30;
opts.Defectmethod = 'PS';
opts.Quadmethod = 'G';
opts.NType = 'LGL';

% if input arguments are provided
% Cart(p,p.nt,opts,opts.Quadmethod,opts.Defectmethod,opts.NType)
if nargin >= 1
    p = varargin{1};
end
if nargin >= 2
    p.nt = varargin{2};
end
if nargin >= 3
    opts = varargin{3};
end
if nargin >= 4
    opts.Quadmethod = varargin{4};
end
if nargin >= 5
    opts.Defectmethod = varargin{5};
end
if nargin >= 6
    opts.NType = varargin{6};
end
if nargin > 6
    warning('too many input arguments...');
end

% set current file name and path
[mpath,mname] = fileparts(mfilename('fullpath'));
opts.mpath = mpath;
opts.mname = mname;

%% tunable parameters
p.t0 = 0; p.tf = 10; % time horizon
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
setup.LB = LB; setup.UB = UB; setup.p = p;

%% solve
[T,U,Y,P,F,p,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = AndersonMoore64_output(T,U,Y,P,F,p,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
AndersonMoore64_plot(T,U,Y,P,F,p,opts,sol)