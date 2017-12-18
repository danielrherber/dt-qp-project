%--------------------------------------------------------------------------
% BrysonHo116.m
% pp. 116-117 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = BrysonHo116(varargin)

% if input arguments are provided
% BrysonHo166(p,p.nt,opts,opts.Quadmethod,opts.Defectmethod,opts.NType)
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
p.tf = 1.5;
p.x0 = 0.1;
p.v0 = 0.5;

%% setup
p.t0 = 0;

% system dynamics
A = [0 1; 0 0]; 
B = [0 0; 1 0]; 

% Lagrange term
L(1).right = 1; % control variables
L(1).matrix = [0,1];

% initial conditions
LB(1).right = 4; % initial states
LB(1).matrix = [p.x0;p.v0];
UB(1).right = 4; % initial states
UB(1).matrix = [p.x0;p.v0];

% final conditions
LB(2).right = 5; % final states
LB(2).matrix = [0;0];
UB(2).right = 5; % final states
UB(2).matrix = [0;0];

% absolute value control bounds
LB(3).right = 1; % control
LB(3).matrix = [-1;-Inf];
UB(3).right = 1; % control
UB(3).matrix = [1;Inf];

% absolute value approximation
Z(1).linear(1).right = 1; % controls
Z(1).linear(1).matrix = [-1;-1];
Z(1).b = 0;
Z(2).linear(1).right = 1; % controls
Z(2).linear(1).matrix = [1;-1];
Z(2).b = 0;

% combine
setup.A = A; setup.B = B; setup.L = L; setup.UB = UB; setup.LB = LB;
setup.Z = Z; setup.p = p;

%% solve
[T,U,Y,P,F,p,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = BrysonHo116_output(T,U,Y,P,F,p,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
BrysonHo116_plot(T,U,Y,P,F,p,opts,sol)