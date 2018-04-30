%--------------------------------------------------------------------------
% BrysonDenham.m
% pp. 120-123 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = BrysonDenham(varargin)

% set p and opts (see BrysonDenham_opts)
% input arguments can be provided in the format 'BrysonDenham(p,opts)'
[p,opts] = DTQP_standardizedinputs(@BrysonDenham_opts,varargin);

%% tunable parameters
p.ell = 1/9;

%% setup
% time horizon
p.t0 = 0; p.tf = 1;

% system dynamics
A = [0 1;0 0]; B = [0;1];

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = 1/2; % 1/2*u^2

% simple bounds
UB(1).right = 4; UB(1).matrix = [0;1]; % initial states
LB(1).right = 4; LB(1).matrix = [0;1];
UB(2).right = 5; UB(2).matrix = [0;-1]; % final states
LB(2).right = 5; LB(2).matrix = [0;-1];
UB(3).right = 2; UB(3).matrix = [p.ell;Inf]; % states

% combine structures
setup.A = A; setup.B = B; setup.L = L; setup.UB = UB; setup.LB = LB; setup.p = p;

%% solve
[T,U,Y,P,F,p,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = BrysonDenham_output(T,U,Y,P,F,p,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
BrysonDenham_plot(T,U,Y,P,F,p,opts,sol)

end
% User options function for BrysonDenham example
function opts = BrysonDenham_opts
% test number
num = 1;

switch num
case 1
    % default parameters
    opts.general.plotflag = 1; % create the plots
    opts.general.saveflag = 0;
    opts.general.displevel = 2;
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 16; % number of nodes
    opts.qp.reorder = 0;
    opts.qp.solver = 'built-in';
    opts.qp.tolerance = 1e-15;
    opts.qp.maxiters = 200;
    opts.qp.disp = 'iter';
case 2
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 4; % number of nodes
end
end