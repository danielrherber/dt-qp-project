%--------------------------------------------------------------------------
% ArbitraryTransfer.m
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = ArbitraryTransfer(varargin)

% set p and opts (see ArbitraryTransfer_opts)
% input arguments can be provided in the format 'ArbitraryTransfer(p,opts)'
[p,opts] = DTQP_standardizedinputs(@ArbitraryTransfer_opts,varargin);

%% tunable parameters
t0 = 0; tf = 10; % time horizon
ny = 18; % number of states
nu = 3; % number of controls

% system dynamics
rng(83233683,'twister') % random number seed
Adensity = rand;
Aeig = -2 + (2 - -2).*rand(ny,1);
p.A = sprandsym(ny,Adensity,Aeig);
p.B = -10 + (10 - -10).*rand(ny,nu);

% initial states
p.x0 = 10*rand(ny,1);

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
LB(1).matrix = p.x0;
UB(1).right = 4; % initial states
UB(1).matrix = p.x0;

% final conditions
LB(2).right = 5; % final states
LB(2).matrix = zeros(ny,1);
UB(2).right = 5; % final states
UB(2).matrix = zeros(ny,1);

% combine
setup.A = A; setup.B = B; setup.L = L;
setup.LB = LB; setup.UB = UB; setup.t0 = t0; setup.tf = tf; setup.p = p;

%% solve
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = ArbitraryTransfer_output(T,U,Y,P,F,in,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
ArbitraryTransfer_plot(T,U,Y,P,F,in,opts,sol)

end
% User options function for ArbitraryTransfer example
function opts = ArbitraryTransfer_opts
% test number
num = 1;

switch num
case 1
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED'; 
    opts.dt.nt = 200; % number of time points
end
end