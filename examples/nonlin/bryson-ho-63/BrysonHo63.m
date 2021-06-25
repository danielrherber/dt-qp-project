%--------------------------------------------------------------------------
% BrysonHo63.m
% p. 63 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = BrysonHo63(varargin)
% input arguments can be provided in the format 'BrysonHo63(p,opts)'

% set local functions
ex_opts = @BrysonHo63_opts; % options function
ex_output = @BrysonHo63_output; % output function
ex_plot = @BrysonHo63_plot; % plot function

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
tf = 3;
p.V = 1;
p.w = 0.8;

%% setup
% time horizon
p.t0 = 0; p.tf = tf;

% number of controls, states, and parameters
n.nu = 1; n.ny = 2;

% system dynamics
element.dynamics = '[V*cos(u1) + w; V*sin(u1)]';
element.parameter_list = 'V w';
element.parameter_values = [p.V p.w];

% Lagrange term
element.lagrange = '-y2*(V*cos(u1) + w)';

% simple bounds
UB(1).right = 4; UB(1).matrix = [0;0]; % initial states
LB(1).right = 4; LB(1).matrix = [0;0];
UB(2).right = 5; UB(2).matrix = [0;0]; % final states
LB(2).right = 5; LB(2).matrix = [0;0];
UB(3).right = 1; UB(3).matrix = 2*pi; % controls
LB(3).right = 1; LB(3).matrix = -2*pi;

% guess
Y0 = [[0,0];[0,0]];
U0 = [[0];[0]];
setup.guess.X = [U0,Y0];

% combine structures
setup.element = element; setup.UB = UB; setup.LB = LB;
setup.t0 = p.t0; setup.tf = p.tf; setup.p = p; setup.n = n;

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
function opts = BrysonHo63_opts
% test number
num = 2;

switch num
    case 1
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100+1; % number of nodes
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 20; % number of nodes
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.dt.meshr.method = 'RICHARDSON-DOUBLING';
    opts.dt.meshr.tolerance = 1e-4;
end

end