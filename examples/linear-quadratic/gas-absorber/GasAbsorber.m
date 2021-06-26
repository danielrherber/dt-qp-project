%--------------------------------------------------------------------------
% GasAbsorber.m
% Section 7.4.4 of R. Luus, Iterative Dynamic Programming. CRC Press, 2000,
% isbn: 1584881488
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = GasAbsorber(varargin)
% input arguments can be provided in the format 'GasAbsorber(auxdata,opts)'

% set local functions
ex_opts = @GasAbsorber_opts; % options function
ex_output = @GasAbsorber_output; % output function
ex_plot = @GasAbsorber_plot; % plot function

% set auxdata and opts (see local_opts)
[auxdata,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
testnum = 3;
switch testnum
    case 1
    n = 10;
    case 2
    n = 25;
    case 3
    n = 50;
    case 4
    n = 100;
    case 5
    n = 200;
end

%% setup
% time horizon
t0 = 0; tf = 15;

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = eye(2); % controls
L(2).left = 2; L(2).right = 2; L(2).matrix = eye(n); % states

% system dynamics
A = spdiags([0.538998*ones(n,1) -1.173113*ones(n,1) 0.634115*ones(n,1)],-1:1,n,n);
B = sparse(n,2);
B(1,1)   = 0.538998;
B(end,2) = 0.634115;

% initial states, simple bounds
y0 = -0.0307-(0.1273-0.0307)/(n-1)*((1:n)'-1);
UB(1).right = 4; UB(1).matrix = y0;
LB(1).right = 4; LB(1).matrix = y0;

% control bounds
LB(2).right = 1; LB(2).matrix = [0 0];

% combine structures
setup.A = A; setup.B = B; setup.L = L; setup.UB = UB; setup.LB = LB;
setup.t0 = t0; setup.tf = tf; setup.auxdata = auxdata;

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
function opts = GasAbsorber_opts
% test number
num = 1;

switch num
case 1
    opts.general.plotflag = 1; % create the plots
    opts.general.displevel = 2;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 300; % number of nodes
    opts.solver.tolerance = 1e-12;
    opts.solver.maxiters = 100;
    opts.solver.display = 'iter';
end

end