%--------------------------------------------------------------------------
% BrysonHo154.m
% pp. 154-155 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = BrysonHo154(varargin)
% input arguments can be provided in the format 'BrysonHo154(auxdata,opts)'

% set local functions
ex_opts = @BrysonHo154_opts; % options function
ex_output = @BrysonHo154_output; % output function
ex_plot = @BrysonHo154_plot; % plot function

% set auxdata and opts (see local_opts)
[auxdata,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
t0 = 0; tf = 100; % time horizon
auxdata.x0 = 100; auxdata.v0 = 1;
auxdata.c1 = 1; auxdata.c2 = 100;

%% setup
% system dynamics
A = [0,1;0,0];
B = [0;1];

% Lagrange term
L(1).left = 1; % control variables
L(1).right = 1; % control variables
L(1).matrix = 1/2; % 1/2*u^2

% Mayer term
M(1).left = 5; % final state variables
M(1).right = 5; % final state variables
M(1).matrix = [auxdata.c2,0;0,auxdata.c1]/2;

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
function opts = BrysonHo154_opts
% test number
num = 1;

switch num
case 1
    % default parameters
    opts = [];
end

end