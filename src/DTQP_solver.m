%--------------------------------------------------------------------------
% DTQP_solver.m
% Obtain the solution to the QP using the selected solver
%--------------------------------------------------------------------------
% NOTE: using quadprog even for linear programs
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [X,F,opts] = DTQP_solver(H,f,A,b,Aeq,beq,lb,ub,opts)

% potentially start the timer
if (opts.displevel > 0) % minimal
    tic % start timer
end

switch opts.solver
    %----------------------------------------------------------------------
    case 'built-in' % built-in MATLAB solvers

        if isempty(H) && isempty(f) % feasibility problem
            % options
            options = optimoptions(@quadprog,'Algorithm','interior-point-convex',...
                'Display',opts.disp,...
                'MaxIterations',opts.maxiters,...
                'ConstraintTolerance',opts.tolerance,...
                'OptimalityTolerance',opts.tolerance,...
                'StepTolerance',opts.tolerance);
            warning('off','optim:quadprog:NullHessian');

            % solve the QP
            [X, F, EXITFLAG] = quadprog([],[],A,b,Aeq,beq,lb,ub,[],options);
            
            warning('on','optim:quadprog:NullHessian');
%         elseif isempty(H) % linear program
%             % options
%             options = optimoptions(@linprog, 'Display', opts.disp,...
%                 'TolFun', 100000*eps, 'algorithm', 'interior-point',...
%                 'MaxIter', 200);
% 
%             % solve the LP
%             [X, F, EXITFLAG] = linprog(f,A,b,Aeq,beq,lb,ub,[],options);

        else % quadratic program
            % options
            options = optimoptions(@quadprog,'Algorithm','interior-point-convex',...
                'Display',opts.disp,...
                'MaxIterations',opts.maxiters,...
                'ConstraintTolerance',opts.tolerance,...
                'OptimalityTolerance',opts.tolerance,...
                'StepTolerance',opts.tolerance);

            % solve the QP
            [X, F, EXITFLAG] = quadprog(H,f,A,b,Aeq,beq,lb,ub,[],options);
            
        end
            
        % return Nan if bad exit flag
        if EXITFLAG < 0
            F = NaN;
        end
    %----------------------------------------------------------------------
end

% end the timer
if (opts.displevel > 0) % minimal
    opts.QPsolvetime = toc; % end the timer
end

% display to the command window
if (opts.displevel > 1) % verbose
    disp(['QP solving time: ', num2str(opts.QPsolvetime), ' s'])
end

end