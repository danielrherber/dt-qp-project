%--------------------------------------------------------------------------
% DTQP_solver_quadprog.m
% Interface to quadprog solver in Matlab Optimization Toolbox
%--------------------------------------------------------------------------
% See https://www.mathworks.com/help/optim/ug/quadprog.html
% NOTE: feasibility/linear programs maybe be solved differently in future
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [X,F,in,opts] = DTQP_solver_quadprog(H,f,A,b,Aeq,beq,lb,ub,in,opts)
% qp options
qp = opts.qp;

% turn off warning
warning('off','optim:quadprog:NullHessian');

if isempty(H) && isempty(f) % feasibility problem (solved using quadprog)
    % options
    options = optimoptions(@quadprog,'Algorithm','interior-point-convex',...
        'Display',qp.disp,...
        'MaxIterations',qp.maxiters,...
        'ConstraintTolerance',qp.tolerance,...
        'OptimalityTolerance',qp.tolerance,...
        'StepTolerance',qp.tolerance);

    % solve the feasbility problem
    [X, F, EXITFLAG] = quadprog([],[],A,b,Aeq,beq,lb,ub,[],options);

elseif isempty(H) % linear program (solved using quadprog)
    % options
    options = optimoptions(@quadprog,'Algorithm','interior-point-convex',...
        'Display',qp.disp,...
        'MaxIterations',qp.maxiters,...
        'ConstraintTolerance',qp.tolerance,...
        'OptimalityTolerance',qp.tolerance,...
        'StepTolerance',qp.tolerance);

    % solve the LP
    [X, F, EXITFLAG] = quadprog([],f,A,b,Aeq,beq,lb,ub,[],options);

else % quadratic program
    % options
    options = optimoptions(@quadprog,'Algorithm','interior-point-convex',...
        'Display',qp.disp,...
        'MaxIterations',qp.maxiters,...
        'ConstraintTolerance',qp.tolerance,...
        'OptimalityTolerance',qp.tolerance,...
        'StepTolerance',qp.tolerance);

    % solve the QP
    [X, F, EXITFLAG] = quadprog(H,f,A,b,Aeq,beq,lb,ub,[],options);

end

% turn on warning
warning('on','optim:quadprog:NullHessian');

% return Nan if bad exit flag
if EXITFLAG < 0
    F = NaN;
end

end