%--------------------------------------------------------------------------
% BrysonHo166.m
% pp. 166-167 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = BrysonHo166(varargin)

% default parameters
opts.plotflag = 1; % create the plots
opts.saveflag = 0;
opts.displevel = 2;
opts.Defectmethod = 'TR';
opts.Quadmethod = 'CTR';
opts.NType = 'ED';
p.nt = 100; % number of nodes
opts.reorder = 0;
opts.solver = 'built-in';
opts.tolerance = 1e-12;
opts.maxiters = 200;
opts.disp = 'iter';

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
p.tf = 20; % time horizon
p.x0 = -0.5; p.v0 = 1; % other

%% setup
p.t0 = 0;
setup.p = p;

% system dynamics
setup.A = [0 1;-1 0]; 
setup.B = [0;1]; 

% Lagrange term
setup.L(1).left = 1; % control variables
setup.L(1).right = 1; % control variables
setup.L(1).matrix(1,1) = 1/2; % 1/2*u.^2

% linear boundary constraints
setup.Y(1).linear(1).right = 4; % initial states
setup.Y(1).linear(1).matrix = [1;0];
setup.Y(1).b = p.x0;
setup.Y(2).linear(1).right = 4; % initial states
setup.Y(2).linear(1).matrix = [0;1];
setup.Y(2).b = p.v0;
setup.Y(3).linear(1).right = 5; % final states
setup.Y(3).linear(1).matrix = [1;0];
setup.Y(3).b = 0;
setup.Y(4).linear(1).right = 5; % final states
setup.Y(4).linear(1).matrix = [0;1];
setup.Y(4).b = 0;

%% solve
[T,U,Y,P,F,p,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = BrysonHo166_output(T,U,Y,P,F,p,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
BrysonHo166_plot(T,U,Y,P,F,p,opts,sol)