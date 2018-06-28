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

% set p and opts (see LQRInhomogeneous_opts)
% input arguments can be provided in the format 'LQRInhomogeneous(p,opts)'
[p,opts] = DTQP_standardizedinputs(@LQRInhomogeneous_opts,varargin);

%% tunable parameters
p.ns = 20; % number of states
p.nu = 10; % number of controls
t0 = 0; % time horizon
tf = 10; % time horizon
p.x0 = linspace(-5,5,p.ns)'; % initial states
rng(393872382,'twister') % specific random seed
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
setup.UB = UB; setup.LB = LB; setup.t0 = t0; setup.tf = tf; setup.p = p;

%% solve
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = LQRInhomogeneous_output(T,U,Y,P,F,in,opts,setup);
if nargout == 1
	varargout{1} = O;
end

%% plot
LQRInhomogeneous_plot(T,U,Y,P,F,in,opts,sol)

end
% User options function for LQRInhomogeneous example
function opts = LQRInhomogeneous_opts
% test number
num = 1;

switch num
case 1
    opts.dt.quadrature = 'CQHS';
    opts.dt.defects = 'HS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100;
case 2
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 100;
case 3
    opts.dt.quadrature = 'CTR';
    opts.dt.defects = 'TR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100;
end
end