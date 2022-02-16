%--------------------------------------------------------------------------
% DTQP_SOLVER_osqp.m
% Interface to osqp solver
%--------------------------------------------------------------------------
% See the following:
% https://osqp.org/docs/get_started/matlab.html
% https://osqp.org/docs/interfaces/matlab.html
% https://github.com/osqp/osqp-matlab
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [X,F,in,opts] = DTQP_SOLVER_osqp(H,f,A,b,Aeq,beq,lb,ub,in,opts)

% qp options
solver = opts.solver;

% initialize
prob = osqp;

% construct matrices
P_ = H;
q_ = f;
A_ = [eye(length(lb));Aeq;A];
l_ = [lb;beq;-inf(size(b))];
u_ = [ub;beq;b];

% settings (default value)
settings.alpha = 1.6; % ADMM overrelaxation parameter (1.6)
settings.rho = 0.1; % ADMM rho step (0.1)
settings.sigma = 0.1; % ADMM sigma step (1e-06)
settings.eps_prim_inf = 1e-5; % primal infeasibility tolerance (1e-04)
settings.eps_dual_inf = 1e-5; % dual infeasibility tolerance (1e-04)
settings.eps_rel = solver.tolerance; % relative tolerance (1e-03)
settings.eps_abs = solver.tolerance; % absolute tolerance (1e-03)
settings.max_iter = solver.maxiters; % maximum number of iterations (4000)
settings.verbose = 1; % print output (True)
settings.scaling = 10; % number of scaling iterations (10)
settings.polish = 1; % perform polishing (False)

% construct the problem
prob.setup(P_, q_, A_, l_, u_, settings);

% solve the problem
res = prob.solve();

% extract
X = res.x;
F = res.info.obj_val;

% store output structure and multipliers
in(1).output = res.info;
opts.lambda = res.y;

end