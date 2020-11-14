%--------------------------------------------------------------------------
% Nagurka.m
% Section 6.4 of R. Luus, Iterative Dynamic Programming. CRC Press, 2000,
% isbn: 1584881488
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = Nagurka(varargin)
% input arguments can be provided in the format 'Nagurka(p,opts)'

% set local functions
ex_opts = @Nagurka_opts; % options function
ex_output = @Nagurka_output; % output function
ex_plot = @Nagurka_plot; % plot function

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
testnum = 2;
switch testnum
    case 1
    n = 6;
    case 2
    n = 20;
    case 3
    n = 50;
end

%% setup
% time horizon
t0 = 0; tf = 1;

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = eye(n); % controls
L(2).left = 2; L(2).right = 2; L(2).matrix = eye(n); % states

% Mayer term
M(1).left = 5; M(1).right = 5; M(1).matrix = diag([10 zeros(1,n-1)]); % final states

% system dynamics
vec = (1:n);
A = [sparse(n-1,1), speye(n-1); ...
    sparse(vec.*(-1).^(vec+1))];
B = speye(n);

% initial states, simple bounds
UB(1).right = 4; UB(1).matrix = vec';
LB(1).right = 4; LB(1).matrix = vec';

% combine structures
setup.A = A; setup.B = B; setup.L = L; setup.M = M;
setup.UB = UB; setup.LB = LB; setup.t0 = t0; setup.tf = tf; setup.p = p;

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
function opts = Nagurka_opts
% test number
num = 1;

switch num
case 1
    opts.general.plotflag = 1; % create the plots
    opts.general.displevel = 2;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 200; % number of nodes
    opts.solver.tolerance = 1e-12;
    opts.solver.maxiters = 100;
    opts.solver.display = 'iter';
case 2
    opts.general.plotflag = 1; % create the plots
    opts.general.displevel = 2;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 20; % number of nodes
    opts.solver.tolerance = 1e-12;
    opts.solver.maxiters = 100;
    opts.solver.display = 'iter';
end

end