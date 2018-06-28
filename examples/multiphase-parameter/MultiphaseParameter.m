%--------------------------------------------------------------------------
% MultiphaseParameter.m
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = MultiphaseParameter(varargin)

% set p and opts (see MultiphaseParameter_opts)
% input arguments can be provided in the format 'MultiphaseParameter(p,opts)'
[p,opts] = DTQP_standardizedinputs(@MultiphaseParameter_opts,varargin);

%% tunable parameters
tf1 = 5; % end of phase 1
tf2 = 25; % end of phase 2
tf3 = 30; % end of phase 3

%% setup
rng(20569180,'twister')

% number of states
n = 6;

% initial time
tf0 = 0;

% initial conditions
Y0 = 10*ones(n,1);

% parameter matrix
G = randi([-1 1],n,1);

% linkage
E1 = [-8;zeros(n-1,1)];
E2 = [0;zeros(n-1,1)];

% combine
p.E1 = E1; p.E2 = E2; p.Y0 = Y0; p.t1 = tf1; p.t2 = tf2; p.t3 = tf3; 

% Lagrange term
L(1).right = 2; % states
L(1).left = 2; % states
L(1).matrix = eye(n);

%--- phase 1
% system dynamics
A1 = gallery('jordbloc',n,-1);

% initial conditions
UB(1).right = 4; % initial states
UB(1).matrix = Y0;
LB(1).right = 4; % initial states
LB(1).matrix = Y0;

% combine
setup(1).A = A1; setup(1).G = G; setup(1).L = L;
setup(1).t0 = tf0; setup(1).tf = tf1; setup(1).p = p;
setup(1).UB = UB; setup(1).LB = LB;

%--- phase 2
% system dynamics
A2 = gallery('hanowa',n,-0.05);

% combine
setup(2).A = A2; setup(2).G = G; setup(2).L = L;
setup(2).t0 = tf1; setup(2).tf = tf2; setup(2).p = p;

%--- phase 3
% system dynamics
A3 = -eye(n);

% combine
setup(3).A = A3; setup(3).L = L;
setup(3).t0 = tf2; setup(3).tf = tf3; setup(3).p = p;

%--- phase 1-2 linkage constraints
clear LY
q = eye(n);
for k = 1:n
    LY(k).left.linear.right = 5; % final states
    LY(k).left.linear.matrix = q(:,k);
    LY(k).right.linear.right = 4; % initial states
    LY(k).right.linear.matrix = -q(:,k);
    LY(k).b = E1(k);
end
LY(end+1).left.linear.right = 3; % parameters
LY(end).left.linear.matrix = 1;
LY(end).right.linear.right = 3; % parameters
LY(end).right.linear.matrix = -1;
LY(end).b = 0;
setup(1).LY = LY;

%--- phase 2-3 linkage constraints
clear LY
q = eye(n);
for k = 1:n
    LY(k).left.linear.right = 5; % final states
    LY(k).left.linear.matrix = q(:,k);
    LY(k).right.linear.right = 4; % initial states
    LY(k).right.linear.matrix = -q(:,k);
    LY(k).b = E2(k);
end
setup(2).LY = LY;

%% solve
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = MultiphaseParameter_output(T,U,Y,P,F,in,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
MultiphaseParameter_plot(T,U,Y,P,F,in,opts,sol)

end
% User options function for MultiphaseParameter example
function opts = MultiphaseParameter_opts
% test number
num = 1;

switch num
case 1
    N = 5;
    opts.dt(1).nt = 5*N;
    opts.dt(2).nt = 20*N;
    opts.dt(3).nt = 5*N;
    opts.dt(1).defects = 'PS';
    opts.dt(1).quadrature = 'G';
    opts.dt(1).mesh = 'LGL';
case 2
    N = 50;
    opts.dt(1).nt = 5*N;
    opts.dt(2).nt = 20*N;
    opts.dt(3).nt = 5*N;
    opts.dt(1).defects = 'TR';
    opts.dt(1).quadrature = 'CTR';
    opts.dt(1).mesh = 'ED';
case 3
    opts.dt.nt = 500;
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
end
end