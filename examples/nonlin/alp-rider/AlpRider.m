%--------------------------------------------------------------------------
% AlpRider.m
% pp. 163-165 in J. T. Betts, "Practical Methods for Optimal Control and
% Estimation Using Nonlinear Programming." Society for Industrial and
% Applied Mathematics, Jan. 2010, doi: 10.1137/1.9780898718577
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan(AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = AlpRider(varargin)
% input arguments can be provided in the format 'AlpRider(p,opts)'

% set local functions
ex_opts = @AlpRider_opts;
ex_output = @AlpRider_output;
ex_plot = @AlpRider_plot;

% set p and opts
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
tf = 20;

%% setup
% time horizon
p.t0 = 0; p.tf = tf;

% number of controls, states, and parameters
n.nu = 2; n.ny = 4;

% system dynamics
str{1} = '[';
str{end+1} = '-10*y1 + u1 + u2; ';
str{end+1} = '-2*y2 + u1 + 2*u2; ';
str{end+1} = '-3*y3 + 5*y4 + u1 - u2; ';
str{end+1} = '5*y3 - 3*y4 + u1 + 3*u2';
str{end+1} = ']';
element.dynamics = horzcat(str{:});

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = diag([1e-2,1e-2]);
L(2).left = 2; L(2).right = 2; L(2).matrix = diag([1e2,1e2,1e2,1e2]);

% inequality constraint
element.g.func = '3*exp(-12*(t-3)^2) + 3*exp(-10*(t-6)^2) + 3*exp(-6*(t-10)^2) + 8*exp(-4*(t-15)^2) + 0.01 - y1^2 - y2^2 - y3^2 - y4^2';
element.g.pathboundary = 1;

% simple bounds
UB(1).right = 4; UB(1).matrix = [2,1,2,1]; % initial states
LB(1).right = 4; LB(1).matrix = [2,1,2,1];
UB(2).right = 5; UB(2).matrix = [2,3,1,-2]; % final states
LB(2).right = 5; LB(2).matrix = [2,3,1,-2];

% guess
Y0 = [[2,1,2,1];[2,3,1,-2]];
U0 = [[0,0];[0,0]];
setup.guess.X = [U0,Y0];

% combine structures
setup.element = element; setup.L = L; setup.UB = UB; setup.LB = LB;
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
function opts = AlpRider_opts
% test number
num = 2;

switch num
case 1
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100; % number of nodes
    opts.solver.maxiters = 300;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.method.olqflag = false;
    opts.dt.meshr.method = 'RICHARDSON-DOUBLING';
    opts.dt.meshr.tolerance = 1e-3;
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 100; % number of nodes
    opts.solver.maxiters = 300;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.method.olqflag = false;
end

end