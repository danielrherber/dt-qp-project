%--------------------------------------------------------------------------
% HDAE.m
% Heat Diffusion Process with Inequality
%--------------------------------------------------------------------------
% pp. 192-195 of J. T. Betts, Practical Methods for Optimal Control and
% Estimation Using Nonlinear Programming. SIAM, 2010, isbn: 9780898716887.
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = HDAE(varargin)
% input arguments can be provided in the format 'HDAE(p,opts)'

% set local functions
ex_opts = @HDAE_opts; % options function
ex_output = @HDAE_output; % output function
ex_plot = @HDAE_plot; % plot function

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
n = 20;
a = 0.5;
b = 0.2;
c = 1;
q1 = 1e-3;
q2 = 1e-3;
tf = 5;

% save to parameter structure
p.n = n;

%% setup
delta = pi/n;

% system dynamics
A = gallery('tridiag',n-1,1,-2,1);
A = A/delta^2;

B = zeros(n-1,2);
B(1,1) = 1;
B(n-1,2) = 1;
B = B/delta^2;

%--- Lagrange term
idx = 0;

% quadratic control
idx = idx + 1;
L(idx).left = 1; % controls;
L(idx).right = 1; % controls;
L(idx).matrix = diag([1/2*delta+q1,1/2*delta+q2]);

% quadratic states
idx = idx + 1;
L(idx).left = 2; % states
L(idx).right = 2; % states
L(idx).matrix = delta*eye(n-1);

%--- simple lower bounds
idx = 0;

% initial states
idx = idx + 1;
LB(idx).right = 4; % initial states
LB(idx).matrix = zeros(n-1,1);

% state inequality constraints
idx = idx + 1;
LB(idx).right = 2; % states
g = cell(1,n-1);
for k = 1:n-1
   g{1,k} = @(t,p)  c*(sin(k*pi/n)*sin(pi*t/5) - a) - b;
end
LB(idx).matrix = g;

% control inequality constraints
idx = idx + 1;
LB(idx).right = 1; % controls
LB(idx).matrix{1,1} = @(t,p)  c*(sin(0*pi/n)*sin(pi*t/5) - a) - b;
LB(idx).matrix{2,1} = @(t,p)  c*(sin(n*pi/n)*sin(pi*t/5) - a) - b;

% initial controls
idx = idx + 1;
LB(idx).right = 1; % controls
LB(idx).matrix{1,1} = @(t) [0;-inf(length(t)-1,1)];
LB(idx).matrix{2,1} = @(t) [0;-inf(length(t)-1,1)];

%--- simple upper bounds
idx = 0;

% initial states
idx = idx + 1;
UB(idx).right = 4; % initial states
UB(idx).matrix = zeros(n-1,1);

% initial controls
idx = idx + 1;
UB(idx).right = 1; % controls
UB(idx).matrix{1,1} = @(t) [0;inf(length(t)-1,1)];
UB(idx).matrix{2,1} = @(t) [0;inf(length(t)-1,1)];

% combine
setup.A = A; setup.B = B; setup.L = L;
setup.LB = LB; setup.UB = UB; setup.tf = tf; setup.p = p;

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
function opts = HDAE_opts
% test number
num = 1;

switch num
case 1
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100;
case 2
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 100; % number of nodes
case 3
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 10;
    opts.solver.tolerance = 1e-15;
    opts.solver.display = 'none';
    opts.dt.meshr.method = 'SS-BETTS';
    opts.dt.meshr.tolerance = 1e-2;
end

end