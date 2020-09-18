%--------------------------------------------------------------------------
% OutputTracking.m
% p. 175 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = OutputTracking(varargin)
% input arguments can be provided in the format 'OutputTracking(p,opts)'

% set local functions
ex_opts = @OutputTracking_opts; % options function
ex_output = @OutputTracking_output; % output function
ex_plot = @OutputTracking_plot; % plot function

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
tf = 50; % final time

ny = 8; % number of states
nu = 3; % number of controls
no = 2; % number of outputs

rng(900911075,'twister') % specific random seed
x0 = ones(ny,1); % initial states
A = sprand(ny,ny,0.7,1); % state matrix
B = sprand(ny,nu,1,1); % input matrix
C = sprand(no,ny,1,1); % output matrix
R = 1e-2*eye(nu); % control penalty
Q = eye(no); % state penalty

% output to track
[o,W] = OutputTracking_o(no);

% save to parameter structure
p.x0 = x0; p.A = A; p.B = B; p.C = C; p.R = R; p.Q = Q; p.W = W; p.o = o;

%% setup
% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = R; % u'*R*u
L(2).left = 2; L(2).right = 2; L(2).matrix = C'*Q*C; % (C*x)'*Q*(C*x)
L(3).left = 0; L(3).right = 2; L(3).matrix = {'prod',o',-2*Q,C}; % -2*o'*Q*C
L(4).left = 0; L(4).right = 0; L(4).matrix = {'prod',o',Q,o}; % o'*Q*o

% initial states, simple bounds
UB(1).right = 4; UB(1).matrix = x0; % initial states
LB(1).right = 4; LB(1).matrix = x0; % initial states

% combine structures
setup.A = A; setup.B = B; setup.L = L;
setup.UB = UB; setup.LB = LB; setup.tf = tf; setup.p = p;

%% solve
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = ex_output(T,U,Y,P,F,in,setup,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
ex_plot(T,U,Y,P,F,in,opts,sol)

end
% User options function for this example
function opts = OutputTracking_opts
% test number
num = 1;

switch num
case 1
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 100; % number of nodes
case 2
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 1000; % number of nodes
case 3
    opts = [];
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 50;
    opts.solver.tolerance = 1e-15;
    opts.solver.display = 'none';
    opts.dt.meshr.method = 'SS-BETTS';
    opts.dt.meshr.tolerance = 1e-3;
end

end