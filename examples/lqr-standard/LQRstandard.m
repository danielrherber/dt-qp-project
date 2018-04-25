%--------------------------------------------------------------------------
% LQRstandard.m
% pp. 148-152 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = LQRstandard(varargin)

% set p and opts (see LQRstandard_opts.m)
% input arguments can be provided in the format 'LQRstandard(p,opts)'
[p,opts] = DTQP_standardizedinputs('LQRstandard_opts',varargin);

%% tunable parameters
p.ns = 20; % number of states
p.nu = 10; % number of controls
p.t0 = 0; % time horizon
p.tf = 10; % time horizon
p.x0 = linspace(-5,5,p.ns)'; % initial states
rng(393872382,'twister') % specific random seed
p.A = sprand(p.ns,p.ns,0.5,1);
p.B = sprand(p.ns,p.nu,1,1);
p.R = eye(p.nu);
p.Q = sprand(p.ns,p.ns,0.2);
p.Q = ((p.Q)*((p.Q)'))/100;
p.M = 10*eye(p.ns); % objective

%% setup
p.t0 = 0;
setup.p = p;

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = p.R; % u'*R*u
L(2).left = 2; L(2).right = 2; L(2).matrix = p.Q; % x'*Q*x

% Mayer term
M(1).left = 5; M(1).right = 5; M(1).matrix = p.M; %xf'*M*xf

% initial states, simple bounds
UB(1).right = 4; UB(1).matrix = p.x0; % initial states
LB(1).right = 4; LB(1).matrix = p.x0;

% combine structures
setup.A = p.A; setup.B = p.B; setup.L = L; setup.M = M;
setup.UB = UB; setup.LB = LB; setup.p = p;

%% solve
[T,U,Y,P,F,p,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = LQRstandard_output(T,U,Y,P,F,p,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
LQRstandard_plot(T,U,Y,P,F,p,opts,sol)