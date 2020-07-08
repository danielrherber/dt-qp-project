%--------------------------------------------------------------------------
% MinimumEnergyTransfer.m
% pp. 21-23 of J.  Klamka, Controllability and Minimum Energy Control,
% Springer, 2019, doi: 10.1007/978-3-319-92540-0
% pp. 163-164 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = MinimumEnergyTransfer(varargin)
% input arguments can be provided in the format 'MinimumEnergyTransfer(p,opts)'

% set local functions
ex_opts = @MinimumEnergyTransfer_opts; % options function
ex_output = @MinimumEnergyTransfer_output; % output function
ex_plot = @MinimumEnergyTransfer_plot; % plot function

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
t0 = 0; tf = 2; % time horizon
ny = 18; % number of states
nu = 6; % number of controls

% system dynamics
rng(83233683,'twister') % random number seed
Adensity = rand;
Aeig = -2 + (2 - -2).*rand(ny,1);
p.A = sprandsym(ny,Adensity,Aeig);
p.B = -10 + (10 - -10).*rand(ny,nu);

% initial states
p.y0 = 10*rand(ny,1);
p.yf = zeros(ny,1);

% check controllability
try
	Co = ctrb(p.A,p.B);
    if rank(Co) ~= ny
        warning('system is not controllable')
    end
catch
    warning('unable to check controllability')
end

%% setup
% system dynamics
A = p.A;
B = p.B;

% Lagrange term
L(1).left = 1; % control variables
L(1).right = 1; % control variables
L(1).matrix = eye(nu);

% initial conditions
LB(1).right = 4; % initial states
LB(1).matrix = p.y0;
UB(1).right = 4; % initial states
UB(1).matrix = p.y0;

% final conditions
LB(2).right = 5; % final states
LB(2).matrix = p.yf;
UB(2).right = 5; % final states
UB(2).matrix = p.yf;

% combine
setup.A = A; setup.B = B; setup.L = L;
setup.LB = LB; setup.UB = UB; setup.t0 = t0; setup.tf = tf; setup.p = p;

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
function opts = MinimumEnergyTransfer_opts
% test number
num = 3;

switch num
case 1
    opts = [];
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100; % number of nodes
case 2
    opts = [];
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 50; % number of nodes
case 3
    opts = [];
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 4; % number of nodes
    opts.solver.display = 'none';
    opts.dt.meshr.method = 'SS-BETTS';
    opts.dt.meshr.tolerance = 1e-8;
end

end