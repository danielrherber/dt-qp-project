%--------------------------------------------------------------------------
% BrysonHo116.m
% pp. 116-117 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = BrysonHo116(varargin)
% input arguments can be provided in the format 'BrysonHo116(auxdata,opts)'

% set local functions
ex_opts = @BrysonHo116_opts; % options function
ex_output = @BrysonHo116_output; % output function
ex_plot = @BrysonHo116_plot; % plot function

% set auxdata and opts (see local_opts)
[auxdata,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
tf = 1.5;
auxdata.x0 = 0.1;
auxdata.v0 = 0.5;

%% setup
t0 = 0;

% system dynamics
A = [0 1; 0 0];
B = [0 0; 1 0];

% Lagrange term
L(1).right = 1; % control variables
L(1).matrix = [0,1];

% initial conditions
LB(1).right = 4; % initial states
LB(1).matrix = [auxdata.x0;auxdata.v0];
UB(1).right = 4; % initial states
UB(1).matrix = [auxdata.x0;auxdata.v0];

% final conditions
LB(2).right = 5; % final states
LB(2).matrix = [0;0];
UB(2).right = 5; % final states
UB(2).matrix = [0;0];

% absolute value control bounds
LB(3).right = 1; % control
LB(3).matrix = [-1;-Inf];
UB(3).right = 1; % control
UB(3).matrix = [1;Inf];

% absolute value approximation
Z(1).linear(1).right = 1; % controls
Z(1).linear(1).matrix = [-1;-1];
Z(1).b = 0;
Z(2).linear(1).right = 1; % controls
Z(2).linear(1).matrix = [1;-1];
Z(2).b = 0;

% combine
setup.A = A; setup.B = B; setup.L = L; setup.UB = UB; setup.LB = LB;
setup.Z = Z; setup.t0 = t0; setup.tf = tf; setup.auxdata = auxdata;

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
function opts = BrysonHo116_opts
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
    opts.solver.display = 'none';
    opts.dt.meshr.method = 'SS-BETTS';
    opts.dt.meshr.tolerance = 1e-7;
end

end