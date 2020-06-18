%--------------------------------------------------------------------------
% DTQP_solver_ipfmincon.m
% Interface to fmincon solver in Matlab Optimization Toolbox
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [X,F,in,opts] = DTQP_solver_ipfmincon(H,f,A,b,Aeq,beq,lb,ub,in,opts)

% extract
X0 = in.X0; dyn = in.dyn; obj = in.obj; cin = in.cin; ceq = in.ceq;

% qp options
qp = opts.qp;

% options cell array
o = {'Display','iter',...
    'MaxIterations',qp.maxiters,...
    'MaxFunctionEvaluations',1e7,...
    'CheckGradients',false,...
    'ConstraintTolerance',qp.tolerance,...
    'OptimalityTolerance',qp.tolerance,...
    'StepTolerance',qp.tolerance};

% algorithm
o{end+1} = 'Algorithm';
o{end+1} = 'interior-point';
% o{end+1} = 'Algorithm';
% o{end+1} = 'sqp';

% use analytic derivatives
if opts.qlin.derivativeflag
    o{end+1} = 'SpecifyConstraintGradient';
    o{end+1} = true;
    o{end+1} = 'SpecifyObjectiveGradient';
    o{end+1} = true;
    o{end+1} = 'HessianFcn';
    o{end+1} = @(x,lambda) DTQP_ipfmincon_hessian(x,lambda,obj,dyn,cin,ceq,H,in,opts);
else
    o{end+1} = 'SpecifyConstraintGradient';
    o{end+1} = false;
    o{end+1} = 'SpecifyObjectiveGradient';
    o{end+1} = false;
    o{end+1} = 'Hessian';
    o{end+1} = 'lbfgs';
end

% TODO: expose these options
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
    [X,F,EXITFLAG,OUTPUT,LAMBDA] = fmincon(@(X) DTQP_ipfmincon_objective(X,obj,in,opts,H,f),...
        X0,A,b,Aeq,beq,lb,ub,[],options);
else
    % solve the NLDO problem
    [X,F,EXITFLAG,OUTPUT,LAMBDA] = fmincon(@(X) DTQP_ipfmincon_objective(X,obj,in,opts,H,f),...
        X0,A,b,Aeq,beq,lb,ub,...
        @(X)DTQP_ipfmincon_constraints(X,dyn,cin,ceq,in,opts),options);
end

% store output structure and multipliers
in(1).output = OUTPUT;
opts.lambda = LAMBDA;

% return nan if negative exit flag
if EXITFLAG < 0
    F = NaN;
end

end