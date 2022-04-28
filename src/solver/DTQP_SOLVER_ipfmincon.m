%--------------------------------------------------------------------------
% DTQP_SOLVER_ipfmincon.m
% Interface to fmincon solver in Matlab Optimization Toolbox
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [X,F,in,opts] = DTQP_SOLVER_ipfmincon(H,f,A,b,Aeq,beq,lb,ub,in,opts)

% extract
X0 = in.X0; dyn = in.dyn; obj = in.obj; cin = in.cin; ceq = in.ceq;

% qp options
solver = opts.solver;

% options cell array
o = {'Display',solver.display,...
    'MaxIterations',solver.maxiters,...
    'MaxFunctionEvaluations',1e7,...
    'CheckGradients',false,...
    'ConstraintTolerance',solver.tolerance,...
    'OptimalityTolerance',solver.tolerance,...
    'StepTolerance',solver.tolerance};

% NOTE: document and update later
% separate feasibility and optimality tolerances
if isfield(opts.solver,'Otolerance')
    o = {'Display',solver.display,...
        'MaxIterations',solver.maxiters,...
        'MaxFunctionEvaluations',1e7,...
        'CheckGradients',false,...
        'ConstraintTolerance',opts.solver.Ftolerance,...
        'OptimalityTolerance',opts.solver.Otolerance,...
        'StepTolerance',opts.solver.Otolerance};
end

% algorithm
o{end+1} = 'Algorithm';
o{end+1} = 'interior-point';
% o{end+1} = 'Algorithm';
% o{end+1} = 'sqp';

% determine options for the selected derivative method
switch opts.method.derivatives
    %----------------------------------------------------------------------
    case 'internal'
    o{end+1} = 'SpecifyConstraintGradient';
    o{end+1} = false;
    o{end+1} = 'SpecifyObjectiveGradient';
    o{end+1} = false;
    o{end+1} = 'Hessian';
    o{end+1} = 'lbfgs';
    %----------------------------------------------------------------------
    case {'real-forward','real-central','complex','symbolic'}
    o{end+1} = 'SpecifyConstraintGradient';
    o{end+1} = true;
    o{end+1} = 'SpecifyObjectiveGradient';
    o{end+1} = true;
    o{end+1} = 'HessianFcn';
    o{end+1} = @(x,lambda) DTQP_NLP_hessian(x,lambda,obj,dyn,cin,ceq,H,in,opts);
    %----------------------------------------------------------------------
end

% TODO: expose these options
o{end+1} = 'FiniteDifferenceType';
o{end+1} = 'central'; % useful when CheckGradients is true
o{end+1} = 'HonorBounds';
o{end+1} = false;
o{end+1} = 'ScaleProblem';
o{end+1} = 'obj-and-constr';
% o{end+1} = 'InitBarrierParam';
% o{end+1} = 1e3;

% combine options
options = optimoptions(@fmincon,o{:});

% TODO: handle other constraint cases
if isempty(dyn)
    % solve the NLDO problem
    [X,F,EXITFLAG,OUTPUT,LAMBDA] = fmincon(@(X) DTQP_NLP_objective(X,obj,in,opts,H,f),...
        X0,A,b,Aeq,beq,lb,ub,[],options);
else
    % solve the NLDO problem
    [X,F,EXITFLAG,OUTPUT,LAMBDA] = fmincon(@(X) DTQP_NLP_objective(X,obj,in,opts,H,f),...
        X0,A,b,Aeq,beq,lb,ub,...
        @(X)DTQP_NLP_constraints(X,dyn,cin,ceq,in,opts),options);
end

% store output structure and multipliers
in(1).output = OUTPUT;
opts.lambda = LAMBDA;

% return nan if negative exit flag
if EXITFLAG < 0
    F = NaN;
end

end