%--------------------------------------------------------------------------
% LQRInhomogeneous.m
% pp. 175-176 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = LQRInhomogeneous(varargin)

% default parameters
opts.plotflag = 1; % create the plots
opts.saveflag = 0;
opts.displevel = 2;
opts.Defectmethod = 'HS';
opts.Quadmethod = 'CQHS';
opts.NType = 'ED';
p.nt = 100; % number of nodes

% if input arguments are provided
% LQRInhomogeneous(p,p.nt,opts,opts.Quadmethod,opts.Defectmethod,opts.NType)
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
if nargin >= 7
	p = varargin{7};
end
if nargin > 7
    warning('too many input arguments...');
end

% set current file name and path
[mpath,mname] = fileparts(mfilename('fullpath'));
opts.mpath = mpath;
opts.mname = mname;

%% tunable parameters
p.ns = 20; % number of states
p.nu = 10; % number of controls
p.t0 = 0; % time horizon
p.tf = 10; % time horizon
p.x0 = linspace(-5,5,p.ns)'; % initial states
rng(393872382) % specific random seed
p.A = sprand(p.ns,p.ns,0.5,1);
p.B = sprand(p.ns,p.nu,1,1);
p.R = eye(p.nu);
p.Q = sprand(p.ns,p.ns,0.2);
p.Q = ((p.Q)*((p.Q)'))/100;
p.M = 10*eye(p.ns); % objective

p.d = cell(p.ns,1);
p.d{1} = @(t) 10*sin(3*t);
p.d{3} = @(t) 2*sin(2*t);
p.d{4} = @(t) 4*sin(6*t);
p.d{4} = @(t) -3*sin(1*t);
p.d{5} = @(t) -1*sin(2*t);
p.d{6} = @(t) 10*sin(3*t);
p.d{7} = @(t) 2*sin(2*t);
p.d{8} = @(t) 4*sin(6*t);
p.d{9} = @(t) -3*sin(1*t);
p.d{10} = @(t) -1*sin(2*t);

%% setup
p.t0 = 0;
setup.p = p;

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = p.R/2; % u'*R*u
L(2).left = 2; L(2).right = 2; L(2).matrix = p.Q/2; % x'*Q*x

% Mayer term
M(1).left = 5; M(1).right = 5; M(1).matrix = p.M/2; %xf'*M*xf

% initial states, simple bounds
UB(1).right = 4; UB(1).matrix = p.x0; % initial states
LB(1).right = 4; LB(1).matrix = p.x0;

% combine structures
setup.A = p.A; setup.B = p.B; setup.d = p.d; setup.L = L; setup.M = M;
setup.UB = UB; setup.LB = LB; setup.p = p;

%% solve
[T,U,Y,P,F,p,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = LQRInhomogeneous_output(T,U,Y,P,F,p,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
LQRInhomogeneous_plot(T,U,Y,P,F,p,opts,sol)