%--------------------------------------------------------------------------
% TransferMinTime.m
% p. 66-69 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = TransferMinTime(varargin)
% input arguments can be provided in the format 'TransferMinTime(auxdata,opts)'

% set local functions
ex_opts = @TransferMinTime_opts; % options function
ex_output = @TransferMinTime_output; % output function
ex_plot = @TransferMinTime_plot; % plot function

% set auxdata and opts (see local_opts)
[auxdata,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
xT  = 0.1405;
m0 = 1;
dm = 0.1405*0.533;
xmu = 1;
rf = 1.525;

%% setup
% time horizon
auxdata.t0 = 0; auxdata.tf = 1;

% number of controls, states, and parameters
n.nu = 1; n.ny = 4; n.np = 1;

% system dynamics
str{1} = '[';
str{end+1} = 'p1*(y3); ';
str{end+1} = 'p1*(y4/y1); ';
str{end+1} = 'p1*(y4^2/y1 - xmu/y1^2 + (xT/(m0 - dm*p1*t))*sin(u1)); ';
str{end+1} = 'p1*(-y3*y4/y1 + (xT/(m0 - dm*p1*t))*cos(u1))';
str{end+1} = ']';
element.dynamics = horzcat(str{:});
element.parameter_list = 'xmu xT m0 dm';
element.parameter_values = [xmu xT m0 dm];

% simple bounds
UB(1).right = 4; UB(1).matrix = [1,0,0,1]; % initial states
LB(1).right = 4; LB(1).matrix = [1,0,0,1];
UB(2).right = 5; UB(2).matrix = [rf,inf,0,sqrt(xmu/rf)]; % final states
LB(2).right = 5; LB(2).matrix = [rf,-inf,0,sqrt(xmu/rf)];
UB(3).right = 3; UB(3).matrix = [100]; % parameters
LB(3).right = 3; LB(3).matrix = [0];
UB(4).right = 1; UB(4).matrix = [2*pi]; % controls
LB(4).right = 1; LB(4).matrix = [0];

% Mayer term
M(1).left = 0; M(1).right = 3; M(1).matrix = 1;

% guess
Y0 = [[1,0,0,1];[1.5,pi,0,0.5]];
U0 = [[0];[2*pi]];
P0 = [[3];[3]];
setup.guess.X = [U0,Y0,P0];

% scaling
setup.scaling(1).right = 1; % controls
setup.scaling(1).matrix = [1];
setup.scaling(2).right = 2; % states
setup.scaling(2).matrix = [1.5,2.5,0.3,1];

% combine structures
setup.element = element; setup.UB = UB; setup.LB = LB; setup.M = M;
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
function opts = TransferMinTime_opts
% test number
num = 2;

switch num
case 1 % ipfmincon method
    opts.general.plotflag = 1; % create the plots
    opts.general.displevel = 2;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 40; % number of nodes
    opts.solver.tolerance = 1e-8;
    opts.solver.maxiters = 2000;
    opts.solver.display = 'iter';
    opts.method.form = 'nonlinearprogram';
case 2
    opts.general.plotflag = 1; % create the plots
    opts.general.displevel = 2;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 20; % number of nodes
    opts.solver.tolerance = 1e-8;
    opts.solver.maxiters = 2000;
    opts.solver.display = 'iter';
    opts.method.form = 'nonlinearprogram';
    opts.dt.meshr.method = 'RICHARDSON-DOUBLING';
    opts.dt.meshr.tolerance = 1e-6;
case 3
    opts.general.plotflag = 1; % create the plots
    opts.general.displevel = 2;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 10; % number of nodes
    opts.solver.tolerance = 1e-8;
    opts.solver.display = 'iter';
    opts.method.form = 'nonlinearprogram';
    opts.dt.meshr.method = 'RICHARDSON-DOUBLING';
    opts.dt.meshr.tolerance = 1e-6;
end

end