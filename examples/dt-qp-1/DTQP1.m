%--------------------------------------------------------------------------
% DTQP1.m
% Complex linear-quadratic dynamic optimization problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = DTQP1(varargin)

% default parameters
opts.plotflag = 1; % create the plots
opts.saveflag = 0;
opts.displevel = 2;
opts.Defectmethod = 'HS';
opts.Quadmethod = 'CQHS';
opts.NType = 'ED';
p.nt = 5000; % number of nodes
opts.reorder = 0;
opts.solver = 'built-in';
opts.tolerance = 1e-15;
opts.maxiters = 100;
opts.disp = 'iter';

% if input arguments are provided
% DTQP1(p,p.nt,opts,opts.Quadmethod,opts.Defectmethod,opts.NType)
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

%% setup
% time horizon
p.t0 = 0; p.tf = 1;

% problem parameters
p.g = @(t) sin(2*pi*t) + 0.5;

% system dynamics
A = [-1,2,0,0;3,-4,0,0;1,2,-1,0;1,0,0,0]; B = [1,0;-1,0;0,1/20;0,0]; G = zeros(4,1);

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = eye(2)/10; % u1^2 + u2^2
L(2).left = 1; L(2).right = 2; L(2).matrix = [1,1,0,0;0,0,0,0]; % u1*y1 + u1*y2
L(3).left = 2; L(3).right = 2; L(3).matrix = zeros(4); L(3).matrix(2,2) = 5; % 5*y2^2
L(4).left = 0; L(4).right = 2; L(4).matrix = {0,@(t) -5*2*p.g(t),0,0}; % ?5*2*g*y2
L(5).left = 0; L(5).right = 0; L(5).matrix{1} = @(t) 5*(p.g(t)).^2; % 5*g^2

% Mayer term
M(1).left = 0; M(1).right = 3; M(1).matrix = 1; % p1
% y2(t0)?y2(tf) = 0, equality constraint
Y(1).linear(1).right = 4; Y(1).linear(1).matrix = [0;1;0;0]; % y2(t0)
Y(1).linear(2).right = 5; Y(1).linear(2).matrix = [0;-1;0;0]; % ?y2(tf)
Y(1).b = 0;

% ?y1 + u2/12 < 0, inequality constraint
Z(1).linear(1).right = 2; Z(1).linear(1).matrix = [-1;0;0;0]; % ?y1
Z(1).linear(2).right = 1; Z(1).linear(2).matrix = [0;1/12]; % u2/12
Z(1).b = 0;

% y3 < p, inequality constraint
Z(2).linear(1).right = 2; Z(2).linear(1).matrix = [0;0;1;0]; % y3
Z(2).linear(2).right = 3; Z(2).linear(2).matrix = -1; % ?p1
Z(2).b = 0;

% initial states, simple bounds
UB(1).right = 4; UB(1).matrix = [2;inf;0.5;0];
LB(1).right = 4; LB(1).matrix = [2;-inf;0.5;0];

% final states, simple bounds
UB(2).right = 5; UB(2).matrix = [inf;inf;inf;0];
LB(2).right = 5; LB(2).matrix = -[inf;inf;inf;0];

% abs(u2) < 10, simple bounds
UB(3).right = 1; UB(3).matrix = [inf;10];
LB(3).right = 1; LB(3).matrix = -[inf;10];

% y2 < g(t), simple bounds
UB(4).right = 2; UB(4).matrix= {inf;@(t) p.g(t);inf;inf};

% combine structures
setup.A = A; setup.B = B; setup.G = G; setup.L = L; setup.M = M;
setup.Y = Y; setup.Z = Z; setup.UB = UB; setup.LB = LB; setup.p = p;

%% solve
[T,U,Y,P,F,p,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = DTQP1_output(T,U,Y,P,F,p,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
DTQP1_plot(T,U,Y,P,F,p,opts,sol)