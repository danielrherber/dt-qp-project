%--------------------------------------------------------------------------
% BrysonHo109.m
% pp. 109-110 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = BrysonHo109(varargin)

% default parameters
opts.plotflag = 1; % create the plots
opts.saveflag = 0;
opts.displevel = 2;
opts.Defectmethod = 'HS';
opts.Quadmethod = 'CQHS';
opts.NType = 'ED';
p.nt = 200; % number of nodes
opts.reorder = 0;
opts.solver = 'built-in';
opts.tolerance = 1e-15;
opts.maxiters = 200;
opts.disp = 'iter';

% if input arguments are provided
% BrysonDenham(p,p.nt,opts,opts.Quadmethod,opts.Defectmethod,opts.NType)
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
p.x0 = 1; p.a = 2; % 1
p.tf = 1; % time horizon
p.g = @(t) t.*cos(20*pi*t) - 1/4;

%% setup
% time horizon
p.t0 = 0; 

% system dynamics
A = 0; B{1,1} = p.g;

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = 1/2; % 1/2*u^2

% Mayer term
M(1).left = 5; M(1).right = 5; M(1).matrix = p.a^2/2; % a^2/2*xf^2

% simple bounds
UB(1).right = 4; UB(1).matrix = p.x0; % initial state
LB(1).right = 4; LB(1).matrix = p.x0;
UB(2).right = 1; UB(2).matrix = 1; % control
LB(2).right = 1; LB(2).matrix = -1;

% combine structures
setup.A = A; setup.B = B; setup.L = L; setup.M = M;
setup.UB = UB; setup.LB = LB; setup.p = p;

%% solve
[T,U,Y,P,F,p,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = BrysonHo109_output(T,U,Y,P,F,p,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
BrysonHo109_plot(T,U,Y,P,F,p,opts,sol)