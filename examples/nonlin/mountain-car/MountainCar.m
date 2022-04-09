%--------------------------------------------------------------------------
% MountainCar.m
% A. A. Melnikov, A. Makmal, H. J. Briegel, "Projective Simulation Applied
% to the Grid-World and the Mountain-Car Problem." ArXiv:1405.5459 [Cs],
% May 2014. arXiv.org, http://arxiv.org/abs/1405.5459
% Also see:
% https://openmdao.github.io/dymos/examples/mountain_car/mountain_car.html
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = MountainCar(varargin)
% input arguments can be provided in the format 'MountainCar(auxdata,opts)'

% set local functions
ex_opts = @MountainCar_opts; % options function
ex_output = @MountainCar_output; % output function
ex_plot = @MountainCar_plot; % plot function

% set auxdata and opts (see local_opts)
[auxdata,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters

%% setup
% time horizon
auxdata.t0 = 0; auxdata.tf = 1;

% number of controls, states, and parameters
n.nu = 1; n.ny = 2; n.np = 1;

% system dynamics
strD{1} = '[';
strD{end+1} = 'p1*y2;';
strD{end+1} = 'p1*(0.001*u1 - 0.0025*cos(3*y1));';
strD{end+1} = ']';
element.dynamics = horzcat(strD{:});
% element.parameter_list = '';
% element.parameter_values = [];

% Mayer term
M(1).right = 3; % parameters
M(1).left = 0; % singleton
M(1).matrix = 1;

% simple bounds
UB(1).right = 4; UB(1).matrix = [-0.5,0]'; % initial states
LB(1).right = 4; LB(1).matrix = [-0.5,0]';
UB(2).right = 1; UB(2).matrix = 1; % controls
LB(2).right = 1; LB(2).matrix = -1;
UB(3).right = 5; UB(3).matrix = [0.5,inf]'; % final states
LB(3).right = 5; LB(3).matrix = [0.5,0]';
UB(4).right = 3; UB(4).matrix = 10000; % parameters
LB(4).right = 3; LB(4).matrix = 0.05;
UB(5).right = 2; UB(5).matrix = [0.5,0.07]'; % states
LB(5).right = 2; LB(5).matrix = [-1.2,-0.07]';

% guess
Y0 = [[-0.5,0];[0.5,1]];
U0 = [[0];[1]];
P0 = [[500];[500]];
setup.guess.X = [U0,Y0,P0];

% scaling
setup.scaling(1).right = 1; % controls
setup.scaling(1).matrix = 1;
setup.scaling(2).right = 2; % states
setup.scaling(2).matrix = [1.2, 0.07];
setup.scaling(3).right = 1; % parameters
setup.scaling(3).matrix = 100;

% combine structures
setup.element = element; setup.M = M; setup.UB = UB; setup.LB = LB;
setup.t0 = auxdata.t0; setup.tf = auxdata.tf; setup.auxdata = auxdata; setup.n = n;

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
function opts = MountainCar_opts
% test number
num = 1;

switch num
case 1
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 1000; % number of nodes
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.solver.maxiters = 20000;
    opts.method.derivatives = 'symbolic';

end

end