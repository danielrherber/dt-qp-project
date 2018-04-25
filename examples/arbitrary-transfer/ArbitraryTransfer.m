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

% set p and opts (see ArbitraryTransfer_opts.m)
% input arguments can be provided in the format 'ArbitraryTransfer(p,opts)'
[p,opts] = DTQP_standardizedinputs('ArbitraryTransfer_opts',varargin);

%% tunable parameters
p.t0 = 0; p.tf = 10; % time horizon
p.ns = 18; % number of states
p.nu = 3; % number of controls

% system dynamics
rng(83233683,'twister') % random number seed
Adensity = rand;
Aeig = -2 + (2 - -2).*rand(p.ns,1);
p.A = sprandsym(p.ns,Adensity,Aeig);
p.B = -10 + (10 - -10).*rand(p.ns,p.nu);

% initial states
p.x0 = 10*rand(p.ns,1);

% check controllability
try
	Co = ctrb(p.A,p.B);
    if rank(Co) ~= p.ns
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
L(1).matrix = eye(p.nu);

% initial conditions
LB(1).right = 4; % initial states
LB(1).matrix = p.x0;
UB(1).right = 4; % initial states
UB(1).matrix = p.x0;

% final conditions
LB(2).right = 5; % final states
LB(2).matrix = zeros(p.ns,1);
UB(2).right = 5; % final states
UB(2).matrix = zeros(p.ns,1);

% combine
setup.A = A; setup.B = B; setup.L = L;
setup.LB = LB; setup.UB = UB; setup.p = p;

%% solve
[T,U,Y,P,F,p,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = ArbitraryTransfer_output(T,U,Y,P,F,p,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
ArbitraryTransfer_plot(T,U,Y,P,F,p,opts,sol)