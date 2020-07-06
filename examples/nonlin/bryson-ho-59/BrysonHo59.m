%--------------------------------------------------------------------------
% BrysonHo59.m
% pp. 59-63 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = BrysonHo59(varargin)
% input arguments can be provided in the format 'BrysonHo59(p,opts)'

% set local functions
ex_opts = @BrysonHo59_opts; % options function
ex_output = @BrysonHo59_output; % output function
ex_plot = @BrysonHo59_plot; % plot function

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
tf = 2.05;
p.a = 1;
p.h = 1;

%% setup
% time horizon
p.t0 = 0; p.tf = tf;

% number of controls, states, and parameters
n.nu = 1; n.ny = 4;

% Mayer term
M(1).left = 0; M(1).right = 5; M(1).matrix = [-1,0,0,0];

% system dynamics
symb.D = '[a*cos(u1); a*sin(u1); y1; y2]';
symb.paramstr = 'a';
symb.param = [p.a];

% simple bounds
UB(1).right = 4; UB(1).matrix = [0,0,0,0]; % initial states
LB(1).right = 4; LB(1).matrix = [0,0,0,0];
UB(2).right = 5; UB(2).matrix = [inf,0,inf,p.h]; % final states
LB(2).right = 5; LB(2).matrix = [-inf,0,-inf,p.h];
UB(3).right = 1; UB(3).matrix = pi; % controls
LB(3).right = 1; LB(3).matrix = -pi;
UB(4).right = 2; UB(4).matrix = [inf,inf,inf,inf]; % states
LB(4).right = 2; LB(4).matrix = [-inf,-inf,-inf,0];

% guess
Y0 = [[0,0,0,0];[0,0,0,p.h]];
U0 = [[0];[0]];
p.guess = [U0,Y0];

% combine structures
setup.symb = symb; setup.M = M; setup.UB = UB; setup.LB = LB;
setup.t0 = p.t0; setup.tf = p.tf; setup.p = p; setup.n = n;

%% solve
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = ex_output(T,U,Y,P,F,in,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
% disp("paused"); pause % for quasilinearization plots
ex_plot(T,U,Y,P,F,in,opts,sol)

end
% User options function for this example
function opts = BrysonHo59_opts
% test number
num = 1;

switch num
case 1 % ipfmincon method
    opts.general.plotflag = 1; % create the plots
    opts.general.saveflag = false;
    opts.general.displevel = 2;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100; % number of nodes
    opts.solver.tolerance = 1e-8;
    opts.solver.maxiters = 30;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 60; % number of nodes
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
case 3 % qlin method
    opts.general.plotflag = 1; % create the plots
    opts.general.displevel = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100; % number of nodes
    opts.solver.maxiters = 200;
    opts.solver.display = 'none';
    opts.method.form = 'qlin';
    opts.method.trustregionflag = true;
    opts.method.sqpflag = false;
    opts.method.delta = 1e9;
    opts.method.improveguess = true;
end

end