%--------------------------------------------------------------------------
% DTQP_SQP_trustregion.m
% Construct and solve a trust region problem for qlin/sqp
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [beq,lb,ub,opts] = DTQP_SQP_trustregion(A,b,Aeq,beq,lb,ub,opts)

% extract
reduction = opts.reduction;
delta = opts.qlin.delta;
xi = opts.qlin.xi;

% check if the objective function decreased
% if reduction < 1e-1
%     delta = delta*xi; % contract
% else
%     delta = delta/xi; % expand
% end

% assign
opts.qlin.delta = delta;

% compute reduced trust bound
xidelta = xi*delta;

% find when lb=ub and store values
Ikeep = (ub == lb);
Vkeep = ub(Ikeep);

% determine the appropriate bounds based on allowed step size
ubt = min(ub,xidelta);
ubt = max(ubt,-xidelta);
lbt = max(lb,-xidelta);
lbt = min(lbt,xidelta);

% ensure simple equality constraints are satisfied
lbt(Ikeep) = Vkeep;
ubt(Ikeep) = Vkeep;

% compute the objective terms (Aeq*V-beq)'*(Aeq*V-beq)
H = Aeq'*Aeq; H = 2*H;
f = -2*Aeq'*beq;

% quadprog options
qp = opts.qp;
tol = 1e-14; % change default tolerance
options = optimoptions(@quadprog,'Algorithm','interior-point-convex',...
    'Display',qp.disp,...
    'MaxIterations',qp.maxiters,...
    'ConstraintTolerance',tol,...
    'OptimalityTolerance',tol,...
    'StepTolerance',tol);

% solve least squares problem
[V,FVAL] = quadprog(H,f',A,b,[],[],lbt,ubt,[],options);

% for debugging
% disp(FVAL)

% redefine beq
beq = Aeq*V;

% determine the appropriate bounds based on allowed step size
ub = min(ub,delta);
ub = max(ub,-delta);
lb = max(lb,-delta);
lb = min(lb,delta);

% ensure simple equality constraints are satisfied
lb(Ikeep) = Vkeep;
ub(Ikeep) = Vkeep;

end