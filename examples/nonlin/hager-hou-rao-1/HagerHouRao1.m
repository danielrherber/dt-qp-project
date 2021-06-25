%--------------------------------------------------------------------------
% HagerHouRao1.m
% W. W. Hager, H. Hou, and A. V. Rao, "Convergence Rate for a Gauss
% Collocation Method Applied to Unconstrained Optimal Control," J Optim
% Theory Appl, vol. 169, no. 3, pp. 801–824, Mar. 2016,
% doi: 10.1007/s10957-016-0929-7
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = HagerHouRao1(varargin)
% input arguments can be provided in the format 'HagerHouRao1(p,opts)'

% set local functions
ex_opts = @HagerHouRao1_opts; % options function
ex_output = @HagerHouRao1_output; % output function
ex_plot = @HagerHouRao1_plot; % plot function

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters

%% setup
% time horizon
p.t0 = 0; p.tf = 1;

% number of controls, states, and parameters
n.nu = 2; n.ny = 2;

% objective
element.lagrange = '2*y1^2*y2^2 + 1.25/(y2^2) + u2/y2 + u1^2 + u2^2';

% system dynamics
element.dynamics = '[y1 + u1/y2 + y1*y2*u2; -y2*(0.5 + y2*u2)]';

% simple bounds
UB(1).right = 4; UB(1).matrix = [1,1]; % initial states
LB(1).right = 4; LB(1).matrix = [1,1];

% guess
Y0 = [[1,1];[1,1]];
U0 = [[0,0];[0,0]];
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
% disp("paused"); pause % for quasilinearization plots
ex_plot(T,U,Y,P,F,in,opts,sol)

end
% User options function for this example
function opts = HagerHouRao1_opts
% test number
num = 2;

switch num
case 1
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 50; % number of nodes
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.solver.tolerance = 1e-12;
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 7; % number of nodes
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.solver.tolerance = 1e-12;

end

end