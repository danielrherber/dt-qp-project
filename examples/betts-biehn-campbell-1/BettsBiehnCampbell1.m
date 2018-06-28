%--------------------------------------------------------------------------
% BettsBiehnCampbell1.m
% J. T. Betts, N. Biehn, and S. L. Campbell, Convergence of Nonconvergent 
% IRK Discretizations of Optimal Control Problems with State Inequality 
% Constraints," SIAM Journal on Scientific Computing, vol. 23, no. 6, 
% pp. 1981-2007, 2002. doi: 10.1137/S1064827500383044
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = BettsBiehnCampbell1(varargin)

% set p and opts (see BettsBiehnCampbell1_opts)
% input arguments can be provided in the format 'BettsBiehnCampbell1(p,opts)'
[p,opts] = DTQP_standardizedinputs(@BettsBiehnCampbell1_opts,varargin);

%% setup
t0 = 34/15; tf = 4;

% system dynamics
A = [0 1; 0 0]; 
B = [0;1]; 

% Lagrange term
L(1).left = 1; % control variables
L(1).right = 1; % control variables
L(1).matrix = 1e-3; % 1e-3*u^2
L(2).left = 2; % state variables
L(2).right = 2; % state variables
L(2).matrix = [1 0; 0 0]; % y1^2

% initial values
LB(1).right = 4; % initial states
LB(1).matrix = [302399/50625; 70304/3375];
UB(1).right = 4; % initial states
UB(1).matrix = [302399/50625; 70304/3375];

% simple bound path constraint
LB(2).right = 2; % states
LB(2).matrix = {@(t) 15 - (t-4).^4;-Inf};

% combine
setup.A = A; setup.B = B; setup.L = L; 
setup.LB = LB; setup.UB = UB; setup.t0 = t0; setup.tf = tf; setup.p = p;

%% solve
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = BettsBiehnCampbell1_output(T,U,Y,P,F,in,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
BettsBiehnCampbell1_plot(T,U,Y,P,F,in,opts,sol)

end
% User options function for BettsBiehnCampbell1 example
function opts = BettsBiehnCampbell1_opts
% test number
num = 2;

switch num
case 1
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 11; % number of nodes
case 2
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 11; % number of nodes
end
end