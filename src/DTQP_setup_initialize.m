%--------------------------------------------------------------------------
% 
% 
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function setup = DTQP_setup_initialize

% for debugging
setup.NEWSTRUCT = true;

%--------------------------------------------------------------------------
% general
%--------------------------------------------------------------------------
% auxiliary data
setup.auxdata = [];

% time horizon
setup.t0 = [];
setup.tf = [];

% counts
counts.nx = []; % old ny
counts.nu = [];
counts.np = [];
counts.nv = []; % old nd
setup.counts = counts;

%--------------------------------------------------------------------------
% linear-quadratic DO problem elements
%--------------------------------------------------------------------------
% objective
lq.lagrange = []; % old L
lq.mayer = []; % old M

% dynamic constraints
dynamics.A = []; % old A
dynamics.B = []; % old B
dynamics.Bp = []; % old G
dynamics.Bv = []; % old d
lq.dynamics = dynamics;

% general constraints
lq.equality = []; % old Y
lq.inequality = []; % old Z
lq.ub = []; % old UB
lq.lb = []; % old LB

% linkage constraints
lq.linkage_equality = []; % old LY
lq.linkage_inequality = []; % old LZ

% assign
setup.lq = lq;

%--------------------------------------------------------------------------
% nonlinear problem elements
%--------------------------------------------------------------------------
% objective
nonlin.lagrange = [];
nonlin.mayer = []; % not implemented

% dynamic constraints
nonlin.dynamics = [];

% general constraints
nonlin.equality = [];
nonlin.inequality = [];

% linkage constraints
nonlin.linkage_equality = []; % old LY, not implemented
nonlin.linkage_inequality = []; % old LZ, not implemented

% symbolic parameters
data.symbols = []; % old parameter_list
data.values = []; % parameter_values
nonlin.data = data;

% assign
setup.nonlin = nonlin;

%--------------------------------------------------------------------------
% solving
%--------------------------------------------------------------------------
% scaling
method.scaling = [];

% initial guess
method.guess = [];

% assign
setup.method = method;

end