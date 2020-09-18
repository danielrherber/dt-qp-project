%--------------------------------------------------------------------------
% LQRstandard.m
% pp. 148-152 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = LQRstandard(varargin)
% input arguments can be provided in the format 'LQRstandard(p,opts)'

% set local functions
ex_opts = @LQRstandard_opts; % options function
ex_output = @LQRstandard_output; % output function
ex_plot = @LQRstandard_plot; % plot function

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

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

%% setup
% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = p.R; % u'*R*u
L(2).left = 2; L(2).right = 2; L(2).matrix = p.Q; % x'*Q*x

% Mayer term
M(1).left = 5; M(1).right = 5; M(1).matrix = p.M; %xf'*M*xf

% initial states, simple bounds
UB(1).right = 4; UB(1).matrix = p.x0; % initial states
LB(1).right = 4; LB(1).matrix = p.x0;

% combine structures
setup.A = p.A; setup.B = p.B; setup.L = L; setup.M = M;
setup.UB = UB; setup.LB = LB; setup.t0 = t0; setup.tf = tf; setup.p = p;

%% solve
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = ex_output(T,U,Y,P,F,in,opts,setup);
if nargout == 1
	varargout{1} = O;
end

%% plot
ex_plot(T,U,Y,P,F,in,opts,sol)

end
% User options function for this example
function opts = LQRstandard_opts
% test number
num = 1;

switch num
case 1
    % default parameters
    opts.general.plotflag = 1; % create the plots
    opts.general.saveflag = 0;
    opts.general.displevel = 2;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 50; % number of nodes
    opts.method.reordervariables = 0;
    opts.solver.function = 'built-in';
    opts.solver.tolerance = 1e-15;
    opts.solver.maxiters = 200;
    opts.solver.display = 'iter';
case 2
    opts.dt.defects = 'ZO';
    opts.dt.quadrature = 'CEF';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 20;
end

end