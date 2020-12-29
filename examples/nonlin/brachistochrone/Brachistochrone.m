%--------------------------------------------------------------------------
% Brachistochrone.m
% pp. 133-134 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = Brachistochrone(varargin)
% input arguments can be provided in the format 'Brachistochrone(p,opts)'

% set local functions
ex_opts = @Brachistochrone_opts; % options function
ex_output = @Brachistochrone_output; % output function
ex_plot = @Brachistochrone_plot; % plot function

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
p.casenum = 2; % see below

switch p.casenum
    case 1 % final x and y specified
    p.g = 10;
    p.xf = 2;
    p.yf = 2;
    case 2 % final x and y specified
    p.g = 10;
    p.xf = 10;
    p.yf = 2;
    case 3 % final x specified
    p.g = 32.1740;
    p.xf = 1;
    case 4 % final x specified with path constraint
    p.g = 32.1740;
    p.xf = 1;
    p.h = 0.1;
    p.theta = atan(0.5);
end

%% setup
% time horizon
p.t0 = 0; p.tf = 1;

% number of controls, states, and parameters
n.nu = 1; n.ny = 3; n.np = 1;

% Mayer term
M(1).left = 0; M(1).right = 3; M(1).matrix = 1; % parameters

% simple bounds
UB(1).right = 4; UB(1).matrix = [0,0,0]; % initial states
LB(1).right = 4; LB(1).matrix = [0,0,0];
UB(3).right = 1; UB(3).matrix = pi; % controls
LB(3).right = 1; LB(3).matrix = -pi;
UB(4).right = 3; UB(4).matrix = 200; % parameters
LB(4).right = 3; LB(4).matrix = 0;

switch p.casenum
    case {1,2}
    % system dynamics
    symb.D = '[p1*y3*sin(u1); p1*y3*cos(u1); p1*g*cos(u1)]';
    symb.paramstr = 'g';
    symb.param = [p.g];

    % simple bounds
    UB(2).right = 5; UB(2).matrix = [p.xf,p.yf,inf]; % final states
    LB(2).right = 5; LB(2).matrix = [p.xf,p.yf,-inf];

    % guess
    Y0 = [[0,0,0];[p.xf,p.yf,0]];
    U0 = [[0];[0]];
    P0 = [[1];[1]];
    setup.guess.X = [U0,Y0,P0];

    case {3,4}
    % system dynamics
    symb.D = '[p1*y3*cos(u1); p1*y3*sin(u1); p1*g*sin(u1)]';
    symb.paramstr = 'g';
    symb.param = [p.g];

    % simple bounds
    UB(2).right = 5; UB(2).matrix = [p.xf,inf,inf]; % final states
    LB(2).right = 5; LB(2).matrix = [p.xf,-inf,-inf];

    % guess
    Y0 = [[0,0,0];[p.xf,0,0]];
    U0 = [[0];[0]];
    P0 = [[1];[1]];
    setup.guess.X = [U0,Y0,P0];

end

if p.casenum == 4
    % state inequality constraint
    symb.cin = 'y2 - y1/2 - 0.1'; % hard-coded at the moment
end

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
function opts = Brachistochrone_opts
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
    opts.dt.nt = 50; % number of nodes
    opts.solver.tolerance = 1e-11;
    opts.solver.maxiters = 1000;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 50; % number of nodes
    opts.solver.tolerance = 1e-11;
    opts.solver.maxiters = 100;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
end

end