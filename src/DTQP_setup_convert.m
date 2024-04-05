%--------------------------------------------------------------------------
% 
% 
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function old = DTQP_setup_convert(new)

if ~isfield(new,'NEWSTRUCT')
    disp("setup conversion not needed")
    old = new;
    return
end

% create temporary variable

%--------------------------------------------------------------------------
% general
%--------------------------------------------------------------------------
% auxiliary data
old.auxdata = new.auxdata;

% time horizon
old.t0 = new.t0;
old.tf = new.tf;

% counts
old.n.ny = new.counts.nx;
old.n.nu = new.counts.nu;
old.n.np = new.counts.np;
old.n.nd = new.counts.nv;

%--------------------------------------------------------------------------
% linear-quadratic DO problem elements
%--------------------------------------------------------------------------
% objective
old.L = new.lq.lagrange;
old.M = new.lq.mayer;

% dynamic constraints
old.A = new.lq.dynamics.A;
old.B = new.lq.dynamics.B;
old.G = new.lq.dynamics.Bp;
old.d = new.lq.dynamics.Bv;

% general constraints
old.Y = new.lq.equality;
old.Z = new.lq.inequality;
old.UB = new.lq.ub;
old.LB = new.lq.lb;

% linkage constraints
old.LY = new.lq.linkage_equality;
old.LZ = new.lq.linkage_inequality;

%--------------------------------------------------------------------------
% nonlinear problem elements
%--------------------------------------------------------------------------
% objective
old.element.lagrange = new.nonlin.lagrange;
old.element.mayer = new.nonlin.mayer; % not implemented

% dynamic constraints
old.element.dynamics = new.nonlin.dynamics;

% general constraints
old.element.h = new.nonlin.equality;
old.element.g = new.nonlin.inequality;

% linkage constraints
old.element.LY = new.nonlin.linkage_equality; % not implemented
old.element.LZ = new.nonlin.linkage_inequality; % not implemented

% symbolic parameters
old.element.parameter_list = new.nonlin.data.symbols;
old.element.parameter_values = new.nonlin.data.values;

%--------------------------------------------------------------------------
% solving
%--------------------------------------------------------------------------
% scaling
old.scaling = new.method.scaling;

% initial guess
old.guess = new.method.guess;

end