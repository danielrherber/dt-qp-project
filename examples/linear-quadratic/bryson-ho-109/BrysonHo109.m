%--------------------------------------------------------------------------
% BrysonHo109.m
% pp. 109-110 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = BrysonHo109(varargin)
% input arguments can be provided in the format 'BrysonHo109(auxdata,opts)'

% set local functions
ex_opts = @BrysonHo109_opts; % options function
ex_output = @BrysonHo109_output; % output function
ex_plot = @BrysonHo109_plot; % plot function

% set auxdata and opts (see local_opts)
[auxdata,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
auxdata.x0 = 1; auxdata.a = 2; % 1
tf = 1; % time horizon

%% setup
% time horizon
t0 = 0;

% system dynamics
A = 0; B{1,1} = @BrysonHo109_g;

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = 1/2; % 1/2*u^2

% Mayer term
M(1).left = 5; M(1).right = 5; M(1).matrix = auxdata.a^2/2; % a^2/2*xf^2

% simple bounds
UB(1).right = 4; UB(1).matrix = auxdata.x0; % initial state
LB(1).right = 4; LB(1).matrix = auxdata.x0;
UB(2).right = 1; UB(2).matrix = 1; % control
LB(2).right = 1; LB(2).matrix = -1;

% combine structures
setup.A = A; setup.B = B; setup.L = L; setup.M = M;
setup.UB = UB; setup.LB = LB; setup.t0 = t0; setup.tf = tf; setup.auxdata = auxdata;

%% solve
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);


%% output
[O,sol] = ex_output(T,U,Y,P,F,in,setup,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
ex_plot(T,U,Y,P,F,in,opts,sol)

end
% User options function for this example
function opts = BrysonHo109_opts
% test number
num = 2;

switch num
case 1
    opts.general.plotflag = 1; % create the plots
    opts.general.saveflag = 0;
    opts.general.displevel = 2;
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100; % number of nodes
    opts.method.reordervariables = 0;
    opts.solver.function = 'built-in';
    opts.solver.tolerance = 1e-15;
    opts.solver.maxiters = 200;
    opts.solver.display = 'iter';
case 2
    opts.general.plotflag = 1; % create the plots
    opts.general.saveflag = 0;
    opts.general.displevel = 2;
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 10; % number of nodes
    opts.method.reordervariables = 0;
    opts.solver.function = 'built-in';
    opts.solver.tolerance = 1e-15;
    opts.solver.maxiters = 200;
    opts.solver.display = 'none';
    opts.dt.meshr.method = 'SS-BETTS';
    opts.dt.meshr.tolerance = 1e-7;
end

end