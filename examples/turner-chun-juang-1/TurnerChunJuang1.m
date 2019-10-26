%--------------------------------------------------------------------------
% TurnerChunJuang1.m
% J. D. Turner, H. M. Chun, J. N. Juang, "Closed-Form Solutions for a Class
% of Optimal Quadratic Tracking Problems", Journal of Optimization Theory 
% and Applications, vol. 47, no. 4, pp. 465-481, Dec. 1985.
% doi: 10.1007/BF00942192
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = TurnerChunJuang1(varargin)
% input arguments can be provided in the format 'TurnerChunJuang1(p,opts)'

% set local functions
ex_opts = @TurnerChunJuang1_opts; % options function
ex_output = @TurnerChunJuang1_output; % output function
ex_plot = @TurnerChunJuang1_plot; % plot function

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
t0 = 0; tf = 10; % time horizon
p.y0 = -5; % initial state
p.em = 1;
p.um = 1;
p.a = 1;
p.S = 5;
p.eta = 10;

%% setup
% system dynamics
A = -1/p.a; 
B = 1;

% Lagrange term
L(1).left = 1; % control variables
L(1).right = 1; % control variables
L(1).matrix = p.um^-2;
L(2).left = 2; % state variables
L(2).right = 2; % state variables
L(2).matrix = p.em^-2;
L(3).left = 0; % singleton
L(3).right = 2; % state variables
L(3).matrix = {@(t) -2*p.em^-2*p.eta};
L(4).left = 0; % singleton
L(4).right = 0; % singleton
L(4).matrix = {@(t) p.em^-2*p.eta.^2};

% Mayer term
M(1).left = 5; % final states
M(1).right = 5; % final states
M(1).matrix = p.S;
M(2).left = 0; % singleton
M(2).right = 5; % final states
M(2).matrix = -2*p.S*p.eta;
M(3).left = 0; % singleton
M(3).right = 0; % singleton
M(3).matrix = p.S*p.eta^2;

% initial states
LB(1).right = 4; % initial states
LB(1).matrix = p.y0;
UB(1).right = 4; % initial states
UB(1).matrix = p.y0;

% combine structures
setup.A = A; setup.B = B; setup.L = L; setup.M = M;
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
function opts = TurnerChunJuang1_opts
% test number
num = 1;

switch num
case 1
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 20;
case 2
    opts.dt.quadrature = 'CQHS';
    opts.dt.defects = 'HS';
    opts.dt.mesh = 'C';
    opts.dt.nt = 100;
case 3
    opts.dt.quadrature = 'CTR';
    opts.dt.defects = 'TR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100;
end

end