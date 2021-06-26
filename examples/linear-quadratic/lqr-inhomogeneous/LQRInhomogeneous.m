%--------------------------------------------------------------------------
% LQRInhomogeneous.m
% pp. 175-176 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = LQRInhomogeneous(varargin)
% input arguments can be provided in the format 'LQRInhomogeneous(auxdata,opts)'

% set local functions
ex_opts = @LQRInhomogeneous_opts; % options function
ex_output = @LQRInhomogeneous_output; % output function
ex_plot = @LQRInhomogeneous_plot; % plot function

% set auxdata and opts (see local_opts)
[auxdata,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
auxdata.ns = 20; % number of states
auxdata.nu = 10; % number of controls
t0 = 0; % time horizon
tf = 10; % time horizon
auxdata.x0 = linspace(-5,5,auxdata.ns)'; % initial states
rng(393872382,'twister') % specific random seed
auxdata.A = sprand(auxdata.ns,auxdata.ns,0.5,1);
auxdata.B = sprand(auxdata.ns,auxdata.nu,1,1);
auxdata.R = eye(auxdata.nu);
auxdata.Q = sprand(auxdata.ns,auxdata.ns,0.2);
auxdata.Q = ((auxdata.Q)*((auxdata.Q)'))/100;
auxdata.M = 10*eye(auxdata.ns); % objective
auxdata.d = LQRInhomogeneous_d(auxdata.ns); % time-varying disturbances

%% setup
% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = auxdata.R/2; % u'*R*u
L(2).left = 2; L(2).right = 2; L(2).matrix = auxdata.Q/2; % x'*Q*x

% Mayer term
M(1).left = 5; M(1).right = 5; M(1).matrix = auxdata.M/2; %xf'*M*xf

% initial states, simple bounds
UB(1).right = 4; UB(1).matrix = auxdata.x0; % initial states
LB(1).right = 4; LB(1).matrix = auxdata.x0;

% combine structures
setup.A = auxdata.A; setup.B = auxdata.B; setup.d = auxdata.d; setup.L = L; setup.M = M;
setup.UB = UB; setup.LB = LB; setup.t0 = t0; setup.tf = tf; setup.auxdata = auxdata;

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
case 4
    opts.dt.quadrature = 'CQHS';
    opts.dt.defects = 'HS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 10;
    opts.solver.tolerance = 1e-15;
    opts.solver.display = 'none';
    opts.dt.meshr.method = 'SS-BETTS';
    opts.dt.meshr.tolerance = 1e-1;
end

end