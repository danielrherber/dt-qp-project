%--------------------------------------------------------------------------
% BrysonHo153.m
% p. 153 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = BrysonHo153(varargin)

% set p and opts (see BrysonHo153_opts)
% input arguments can be provided in the format 'BrysonHo153(p,opts)'
[p,opts] = DTQP_standardizedinputs(@BrysonHo153_opts,varargin);

%% tunable parameters
p.t0 = 0; p.tf = 1; % time horizon
p.c = 2; p.x0 = 1; % other

%% setup
% system dynamics
A = 0; 
B = 1; 

% Lagrange term
L(1).left = 1; % control variables
L(1).right = 1; % control variables
L(1).matrix = 1/2; % 1/2*u^2

% Mayer term
M(1).right = 5; % final states
M(1).left = 5; % final states
M(1).matrix = p.c/2; % c/2*yf^2

% initial conditions
LB(1).right = 4; % initial states
LB(1).matrix = p.x0;
UB(1).right = 4; % initial states
UB(1).matrix = p.x0;

% combine
setup.A = A; setup.B = B; setup.L = L; setup.M = M;
setup.LB = LB; setup.UB = UB; setup.p = p;

%% solve
[T,U,Y,P,F,p,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = BrysonHo153_output(T,U,Y,P,F,p,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
BrysonHo153_plot(T,U,Y,P,F,p,opts,sol)

end
% User options function for BrysonHo153 example
function opts = BrysonHo153_opts
% test number
num = 1;

switch num
case 1
    % default parameters
    opts = [];
end
end