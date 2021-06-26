%--------------------------------------------------------------------------
% BrysonHo64.m
% p. 64 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = BrysonHo64(varargin)
% input arguments can be provided in the format 'BrysonHo64(auxdata,opts)'

% set local functions
ex_opts = @BrysonHo64_opts; % options function
ex_output = @BrysonHo64_output; % output function
ex_plot = @BrysonHo64_plot; % plot function

% set auxdata and opts (see local_opts)
[auxdata,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
auxdata.a = 1;
auxdata.l = 0.45;

%% setup
% number of controls, states, and parameters
n.nu = 1; n.ny = 1;

% time horizon
auxdata.t0 = -auxdata.l; auxdata.tf = auxdata.l;

% objective function
element.lagrange = '2*pi*y1*sqrt(1+u1^2)';
% element.lagrange = 'y1^2*(1+u1^2)';

% dynamics
setup.A = 0;
setup.B = 1;

% bounds
UB(1).right = 4; UB(1).matrix = auxdata.a;
LB(1).right = 4; LB(1).matrix = auxdata.a;
UB(2).right = 5; UB(2).matrix = auxdata.a;
LB(2).right = 5; LB(2).matrix = auxdata.a;
LB(3).right = 2; LB(3).matrix = 0;

%% setup
setup.element = element; setup.UB = UB; setup.LB = LB;
setup.t0 = auxdata.t0; setup.tf = auxdata.tf; setup.auxdata = auxdata; setup.n = n;

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
function opts = BrysonHo64_opts
% test number
num = 2;

switch num
case 1
    opts.general.plotflag = 1; % create the plots
    opts.general.saveflag = false;
    opts.general.displevel = 2;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 30; % number of nodes
    opts.solver.tolerance = 1e-8;
    opts.solver.maxiters = 200;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 10; % number of nodes
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
end

end