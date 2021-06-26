%--------------------------------------------------------------------------
% BrysonHo156.m
% p. 156 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = BrysonHo156(varargin)
% input arguments can be provided in the format 'BrysonHo156(auxdata,opts)'

% set local functions
ex_opts = @BrysonHo156_opts; % options function
ex_output = @BrysonHo156_output; % output function
ex_plot = @BrysonHo156_plot; % plot function

% set auxdata and opts (see local_opts)
[auxdata,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
t0 = 0; tf = 10; % time horizon
auxdata.x0 = 10; auxdata.v0 = 1; auxdata.c = 1;
auxdata.omega = 1;

%% setup
% system dynamics
A = [0,1;-auxdata.omega^2,0];
B = [0;1];

% Lagrange term
L(1).left = 1; % control variables
L(1).right = 1; % control variables
L(1).matrix = 1/2; % 1/2*u^2

% Mayer term
M(1).left = 5; % final state variables
M(1).right = 5; % final state variables
M(1).matrix = [auxdata.c/2,0;0,0];

% initial conditions
LB(1).right = 4; % initial states
LB(1).matrix = [auxdata.x0;auxdata.v0];
UB(1).right = 4; % initial states
UB(1).matrix = [auxdata.x0;auxdata.v0];

% combine
setup.A = A; setup.B = B; setup.L = L; setup.M = M;
setup.LB = LB; setup.UB = UB; setup.t0 = t0; setup.tf = tf; setup.auxdata = auxdata;

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
function opts = BrysonHo156_opts
% test number
num = 2;

switch num
case 1
    % default parameters
    opts = [];
case 2
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 5;
    opts.solver.tolerance = 1e-15;
    opts.solver.maxiters = 200;
    opts.solver.display = 'none';
    opts.dt.meshr.method = 'SS-BETTS';
    opts.dt.meshr.tolerance = 1e-8;
end

end