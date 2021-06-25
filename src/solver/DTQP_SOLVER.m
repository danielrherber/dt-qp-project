%--------------------------------------------------------------------------
% DTQP_SOLVER.m
% Obtain the solution to the DO problem using the selected solver
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [X,F,in,opts] = DTQP_SOLVER(H,f,A,b,Aeq,beq,lb,ub,in,opts)

% extract
displevel = opts.general.displevel;

% (potentially) stop create timer and start qpsolver timer
if (displevel > 0) % minimal
    ctime = toc(opts.timer.t3);
    opts.timer.create = opts.timer.create + ctime; % add
    opts.timer.t3 = tic; % start timer
end

% (potentially) display to the command window
if (displevel > 1) % verbose
    flag = 'creation-time'; DTQP_commandWindowTasks %#ok<NASGU>
end

switch upper(opts.solver.function)
    %----------------------------------------------------------------------
    case {'QUADPROG','BUILT-IN'} % built-in MATLAB qp solver
    [X,F,in,opts] = DTQP_SOLVER_quadprog(H,f,A,b,Aeq,beq,lb,ub,in,opts);
    %----------------------------------------------------------------------
    case 'CVX' % http://cvxr.com/cvx/
    [X,F,in,opts] = DTQP_SOLVER_cvx(H,f,A,b,Aeq,beq,lb,ub,in,opts);
    %----------------------------------------------------------------------
    case 'QPOASES' % https://github.com/coin-or/qpOASES
    [X,F,in,opts] = DTQP_SOLVER_qpoases(H,f,A,b,Aeq,beq,lb,ub,in,opts);
    %----------------------------------------------------------------------
    case 'IPFMINCON' % MATLAB fmincon solver with interior point algorithm
    [X,F,in,opts] = DTQP_SOLVER_ipfmincon(H,f,A,b,Aeq,beq,lb,ub,in,opts);
    %----------------------------------------------------------------------
    otherwise
    error('DTQP error: solver not found')
end

% (potentially) stop qpsolver timer and start create qpsolver timer
if (displevel > 0) % minimal
    opts.timer.qpsolver = opts.timer.qpsolver + toc(opts.timer.t3); % add
    opts.timer.t3 = tic; % start timer
end

% (potentially) display to the command window
if (displevel > 1) % verbose
    flag = 'solver-time'; DTQP_commandWindowTasks %#ok<NASGU>
end

end