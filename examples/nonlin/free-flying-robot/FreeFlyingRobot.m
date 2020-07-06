%--------------------------------------------------------------------------
% FreeFlyingRobot.m
% 326-330 of J. T. Betts, "Practical Methods for Optimal Control and
% Estimation Using Nonlinear Programming*." Society for Industrial and
% Applied Mathematics, Jan. 2010, doi: 10.1137/1.9780898718577
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = FreeFlyingRobot(varargin)
% input arguments can be provided in the format 'FreeFlyingRobot(p,opts)'

% set local functions
ex_opts = @FreeFlyingRobot_opts;
ex_output = @FreeFlyingRobot_output;
ex_plot = @FreeFlyingRobot_plot;

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
tf = 12;
a = 0.2; b = 0.2;

%% setup
% time horizon
p.t0 = 0; p.tf = tf;

% number of controls, states, and parameters
n.nu = 4; n.ny = 6;

% system dynamics
str{1} = '[';
str{end+1} = 'y4;';
str{end+1} = 'y5;';
str{end+1} = 'y6;';
str{end+1} = '(u1-u2+u3-u4)*cos(y3);';
str{end+1} = '(u1-u2+u3-u4)*sin(y3);';
str{end+1} = 'a*(u1-u2) - b*(u3-u4)';
str{end+1} = ']';
symb.D = horzcat(str{:});
symb.paramstr = 'a b';
symb.param = [a b];

% Lagrange term
% symb.Ob = 'u1+u2+u3+u4';
L(1).left = 0; L(1).right = 1; L(1).matrix = [1 1 1 1];
setup.L = L;

% inequality constraints
symb.cin.func = '[u1+u2-1;u3+u4-1]';
symb.cin.pathboundary = [1 1];

% simple bounds
UB(1).right = 4; UB(1).matrix = [-10;-10;pi/2;0;0;0]; % initial states
LB(1).right = 4; LB(1).matrix = [-10;-10;pi/2;0;0;0];
UB(2).right = 5; UB(2).matrix = [0;0;0;0;0;0]; % final states
LB(2).right = 5; LB(2).matrix = [0;0;0;0;0;0];
UB(3).right = 1; UB(3).matrix = 0.3+[1;1;1;1]; % controls
LB(3).right = 1; LB(3).matrix = [0;0;0;0];

% guess
Y0 = [[-10,-10,pi/2,0,0,0];[0,0,0,0,0,0]];
U0 = [[0,0,0,0];[0,0,0,0]];
p.guess = [U0,Y0];

% combine structures
setup.symb = symb; setup.UB = UB; setup.LB = LB;
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
function opts = FreeFlyingRobot_opts
% test number
num = 1;

switch num
case 1
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 300; % number of nodes
    opts.solver.display = 'iter';
    opts.solver.maxiters = 2000;
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.method.olqflag = true;
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 100; % number of nodes
    opts.solver.display = 'iter';
    opts.solver.maxiters = 2000;
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.method.olqflag = true;
end

end