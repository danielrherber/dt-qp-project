%--------------------------------------------------------------------------
% DTQP_solver.m
% Obtain the solution to the QP using the selected solver
%--------------------------------------------------------------------------
% NOTE: using quadprog even for linear programs
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [X,F,in,opts] = DTQP_solver(H,f,A,b,Aeq,beq,lb,ub,in,opts)

% potentially start the timer
if (opts.general.displevel > 0) % minimal
    tstart = toc; % start timer
end

switch upper(opts.qp.solver)
    %----------------------------------------------------------------------
    case {'QUADPROG','BUILT-IN'} % built-in MATLAB solvers
        [X,F,in,opts] = DTQP_solver_quadprog(H,f,A,b,Aeq,beq,lb,ub,in,opts);
    %----------------------------------------------------------------------
    case 'CVX' % 
        [X,F,in,opts] = DTQP_solver_cvx(H,f,A,b,Aeq,beq,lb,ub,in,opts);
    %----------------------------------------------------------------------
    case 'QPOASES' % 
        [X,F,in,opts] = DTQP_solver_qpoases(H,f,A,b,Aeq,beq,lb,ub,in,opts);
    %----------------------------------------------------------------------
    otherwise
        error('Error: solver not found')
end

% end the timer
if (opts.general.displevel > 0) % minimal
    in(end).QPsolvetime = toc - tstart; % end the timer
end

% display to the command window
if (opts.general.displevel > 1) % verbose
    disp(['QP solving time: ', num2str(in(end).QPsolvetime), ' s'])
end

end