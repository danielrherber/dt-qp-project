%--------------------------------------------------------------------------
% BrysonHo166.m
% pp. 166-167 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = BrysonHo166(varargin)
% input arguments can be provided in the format 'BrysonHo166(auxdata,opts)'

% set local functions
ex_opts = @BrysonHo166_opts; % options function
ex_output = @BrysonHo166_output; % output function
ex_plot = @BrysonHo166_plot; % plot function

% set auxdata and opts (see local_opts)
[auxdata,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
tf = 20; % time horizon
auxdata.x0 = -0.5; auxdata.v0 = 1; % other

%% setup
t0 = 0;

% system dynamics
A = [0 1;-1 0];
B = [0;1];

% Lagrange term
L(1).left = 1; % control variables
L(1).right = 1; % control variables
L(1).matrix = 1/2; % 1/2*u^2

% initial conditions
LB(1).right = 4; % initial states
LB(1).matrix = [auxdata.x0;auxdata.v0];
UB(1).right = 4; % initial states
UB(1).matrix = [auxdata.x0;auxdata.v0];

% final conditions
LB(2).right = 5; % final states
LB(2).matrix = [0;0];
UB(2).right = 5; % final states
UB(2).matrix = [0;0];

% combine
setup.A = A; setup.B = B; setup.L = L;
setup.LB = LB; setup.UB = UB; setup.t0 = t0; setup.tf = tf; setup.auxdata = auxdata;

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
function opts = BrysonHo166_opts
% test number
num = 2;

switch num
case 1
    % default parameters
    opts.general.plotflag = 1; % create the plots
    opts.general.saveflag = 0;
    opts.general.displevel = 2;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100; % number of nodes
    opts.method.reordervariables = 0;
    opts.solver.function = 'built-in';
    opts.solver.tolerance = 1e-12;
    opts.solver.maxiters = 200;
case 2
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 5;
    opts.method.reordervariables = 0;
    opts.solver.function = 'built-in';
    opts.solver.tolerance = 1e-15;
    opts.solver.maxiters = 200;
    opts.solver.display = 'none';
    opts.dt.meshr.method = 'SS-BETTS';
    opts.dt.meshr.tolerance = 1e-4;
case 3
    opts.dt.defects = 'PS-MI';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 10; % number of intervals
    opts.dt.nn = 12; % polynomial order in each interval
end

end