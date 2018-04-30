%--------------------------------------------------------------------------
% BrysonHo166.m
% pp. 166-167 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = BrysonHo166(varargin)

% set p and opts (see BrysonHo166_opts)
% input arguments can be provided in the format 'BrysonHo166(p,opts)'
[p,opts] = DTQP_standardizedinputs(@BrysonHo166_opts,varargin);

%% tunable parameters
p.tf = 20; % time horizon
p.x0 = -0.5; p.v0 = 1; % other

%% setup
p.t0 = 0;

% system dynamics
A = [0 1;-1 0]; 
B = [0;1]; 

% Lagrange term
L(1).left = 1; % control variables
L(1).right = 1; % control variables
L(1).matrix = 1/2; % 1/2*u^2

% initial conditions
LB(1).right = 4; % initial states
LB(1).matrix = [p.x0;p.v0];
UB(1).right = 4; % initial states
UB(1).matrix = [p.x0;p.v0];

% final conditions
LB(2).right = 5; % final states
LB(2).matrix = [0;0];
UB(2).right = 5; % final states
UB(2).matrix = [0;0];

% combine
setup.A = A; setup.B = B; setup.L = L;
setup.LB = LB; setup.UB = UB; setup.p = p;

%% solve
[T,U,Y,P,F,p,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = BrysonHo166_output(T,U,Y,P,F,p,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
BrysonHo166_plot(T,U,Y,P,F,p,opts,sol)

end
% User options function for BrysonHo166 example
function opts = BrysonHo166_opts
% test number
num = 1;

switch num
case 1
    % default parameters
    opts.general.plotflag = 1; % create the plots
    opts.general.saveflag = 0;
    opts.general.displevel = 2;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100; % number of nodes
    opts.qp.reorder = 0;
    opts.qp.solver = 'built-in';
    opts.qp.tolerance = 1e-12;
    opts.qp.maxiters = 200;
    opts.qp.disp = 'iter';
end
end