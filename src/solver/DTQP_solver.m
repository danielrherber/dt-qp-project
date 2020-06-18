%--------------------------------------------------------------------------
% DTQP_solver.m
% Obtain the solution to the DO problem using the selected solver
%--------------------------------------------------------------------------
% NOTE: using quadprog even for linear programs
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [X,F,in,opts] = DTQP_solver(H,f,A,b,Aeq,beq,lb,ub,in,opts)

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
    disp(strcat("-> Creation time: ",string(ctime)," s"))
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
    case 'IPFMINCON' % MATLAB fmincon solver with interior point algorithm
        [X,F,in,opts] = DTQP_solver_ipfmincon(H,f,A,b,Aeq,beq,lb,ub,in,opts);
    %----------------------------------------------------------------------
    otherwise
        error('Error: solver not found')
end

% (potentially) stop qpsolver timer and start create qpsolver timer
if (displevel > 0) % minimal
    opts.timer.qpsolver = opts.timer.qpsolver + toc(opts.timer.t3); % add
    opts.timer.t3 = tic; % start timer
end

% (potentially) display to the command window
if (displevel > 1) % verbose
    disp(['-> QP solver time: ', num2str(opts.timer.qpsolver), ' s'])
end

end