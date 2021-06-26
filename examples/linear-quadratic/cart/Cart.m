%--------------------------------------------------------------------------
% Cart.m
% pp. 124-126 of D. H. Ballard, An Introduction to Natural Computation. 
% MIT Press, 1999, isbn: 9780262522588
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = Cart(varargin)
% input arguments can be provided in the format 'Cart(auxdata,opts)'

% set local functions
ex_opts = @Cart_opts; % options function
ex_output = @Cart_output; % output function
ex_plot = @Cart_plot; % plot function

% set auxdata and opts (see local_opts)
[auxdata,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
t0 = 0;

%% setup
tf = 1; % time horizon

% system dynamics
A = [0 1; 0 -1]; 
B = [0;1];

% Lagrange term
L(1).left = 1; % control variables
L(1).right = 1; % control variables
L(1).matrix(1,1) = 1/2; % 1/2*u.^2

% Mayer term
M(1).left = 0; % singleton
M(1).right = 5; % final states
M(1).matrix = [-1,0];

% initial states
LB(1).right = 4; % initial states
LB(1).matrix = [0;0];
UB(1).right = 4; % initial states
UB(1).matrix = [0;0];

% combine structures
setup.A = A; setup.B = B; setup.L = L; setup.M = M; 
setup.LB = LB; setup.UB = UB; setup.t0 = t0; setup.tf = tf;

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
function opts = Cart_opts
% test number
num = 1;

switch num
case 1
    opts.dt.nt = 100; % number of nodes
end

end