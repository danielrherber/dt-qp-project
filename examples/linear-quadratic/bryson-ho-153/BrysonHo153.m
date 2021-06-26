%--------------------------------------------------------------------------
% BrysonHo153.m
% p. 153 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = BrysonHo153(varargin)
% input arguments can be provided in the format 'BrysonHo153(auxdata,opts)'

% set local functions
ex_opts = @BrysonHo153_opts; % options function
ex_output = @BrysonHo153_output; % output function
ex_plot = @BrysonHo153_plot; % plot function

% set auxdata and opts (see local_opts)
[auxdata,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
t0 = 0; tf = 1; % time horizon
auxdata.c = 2; auxdata.x0 = 1; % other

%% setup
% system dynamics
A = 0;
B = 1;

% Lagrange term
L(1).left = 1; % control variables
L(1).right = 1; % control variables
L(1).matrix = 1/2; % 1/2*u^2

% Mayer term
M(1).right = 5; % final states
M(1).left = 5; % final states
M(1).matrix = auxdata.c/2; % c/2*yf^2

% initial conditions
LB(1).right = 4; % initial states
LB(1).matrix = auxdata.x0;
UB(1).right = 4; % initial states
UB(1).matrix = auxdata.x0;

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
function opts = BrysonHo153_opts
% test number
num = 1;

switch num
case 1
    % default parameters
    opts = [];
end

end