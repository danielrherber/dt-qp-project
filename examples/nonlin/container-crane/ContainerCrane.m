%--------------------------------------------------------------------------
% ContainerCrane.m
% D. Augustin and H. Maurer, "Sensitivity Analysis and Real-Time Control of
% a Container Crane under State Constraints," in Online Optimization of
% Large Scale Systems, Springer Berlin Heidelberg, 2001, pp. 69–82
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan(AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = ContainerCrane(varargin)
% input arguments can be provided in the format 'ContainerCrane(p,opts)'

% set local functions
ex_opts = @ContainerCrane_opts;
ex_output = @ContainerCrane_output;
ex_plot = @ContainerCrane_plot;

% set p and opts
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
tf = 9;

%% setup
% time horizon
p.t0 = 0; p.tf = tf;

% number of controls, states, and parameters
n.nu = 2; n.ny = 6;

% system dynamics
str{1} = '[';
str{end+1} = 'y4; ';
str{end+1} = 'y5; ';
str{end+1} = 'y6; ';
str{end+1} = 'u1+17.2656*y3; ';
str{end+1} = 'u2; ';
str{end+1} = '-(u1 + 27.0756*y3+2*y5*y6)/y2';
str{end+1} = ']';
symb.D = horzcat(str{:});

% Lagrange term
if ~isfield(opts,'qlin') || ~isfield(opts.qlin,'olqflag') || opts.qlin.olqflag
    L(1).left = 1; L(1).right = 1; L(1).matrix = diag([0.005,0.005]);
    L(2).left = 2; L(2).right = 2; L(2).matrix = diag([0,0,0.5,0,0,0.5]);
    setup.L = L;
else
    symb.Ob = '0.5*(y3^2 + y6^2 + 0.01*(u1^2+u2^2))';
end

% simple bounds
UB(1).right = 4; UB(1).matrix = [0;22;0;0;-1;0]; % initial states
LB(1).right = 4; LB(1).matrix = [0;22;0;0;-1;0];
UB(2).right = 5; UB(2).matrix = [10;14;0;2.5;0;0]; % final states
LB(2).right = 5; LB(2).matrix = [10;14;0;2.5;0;0];
UB(3).right = 2; UB(3).matrix = [inf;inf;inf;2.5;1;inf]; % states
LB(3).right = 2; LB(3).matrix = [-inf;-inf;-inf;-2.5;-1;-inf];
UB(4).right = 1; UB(4).matrix = [2.8337,0.71265]; % controls
LB(4).right = 1; LB(4).matrix = [-2.8337,-0.80865];

% guess
Y0 = [[0,22,0,0,-1,0];[10,14,0,2.5,0,0]];
U0 = [[2.8337,0.71265];[-2.8337,-0.80865]];
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
function opts = ContainerCrane_opts
% test number
num = 1;

switch num
case 1
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 200; % number of nodes
    opts.qp.maxiters = 4000;
    opts.qp.disp = 'iter';
    opts.qp.solver = 'ipfmincon';
    opts.qlin.method = 'ipfmincon';
    opts.qlin.olqflag = true;
    opts.qlin.derivativeflag = true;
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 30; % number of nodes
    opts.qp.maxiters = 4000;
    opts.qp.disp = 'iter';
    opts.qp.solver = 'ipfmincon';
    opts.qlin.method = 'ipfmincon';
    opts.qlin.olqflag = true;
    opts.qlin.derivativeflag = true;
case 3
    opts.general.displevel = 1;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 200; % number of nodes
    opts.qlin.method = 'qlin';
    opts.qp.disp = 'none';
    opts.qlin.trustregionflag = false;
    opts.qlin.improveX0flag = false;
end

end