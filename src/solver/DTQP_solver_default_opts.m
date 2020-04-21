%--------------------------------------------------------------------------
% DTQP_solver_default_opts.m
% Defaults options for the QP solvers
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function opts = DTQP_solver_default_opts(opts)

switch upper(opts.qp.solver)
    %----------------------------------------------------------------------
    case {'QUADPROG','BUILT-IN'} % see DTQP_solver_quadprog.m
    % tolerance
    if ~isfield(opts.qp,'tolerance')
        opts.qp.tolerance = 1e-12;
    end

    % maximum iterations
    if ~isfield(opts.qp,'maxiters')
        opts.qp.maxiters = 200;
    end

    % display level in the optimization routine
    if ~isfield(opts.qp,'disp')
        opts.qp.disp = 'none'; % none
        if opts.general.displevel > 1 % verbose
            opts.qp.disp = 'iter'; % iterations
        end
        if ~strcmpi(opts.dt(1).meshr.method,'none')
            opts.qp.disp = 'none'; % iterations
        end
    end
    %----------------------------------------------------------------------
    case 'CVX' % see DTQP_solver_cvx.m
    %----------------------------------------------------------------------
    case 'QPOASES' % see DTQP_solver_qpoases.m
    % maximum iterations
    if ~isfield(opts.qp,'maxIter')
        opts.qp.maxIter = -1; % value chosen heuristically
    end

    % maximum CPU time in seconds
    if ~isfield(opts.qp,'maxCpuTime')
        opts.qp.maxCpuTime = -1; % only iteration limit is used
    end

    % display level in the optimization routine
    if ~isfield(opts.qp,'printLevel')
        if opts.general.displevel > 1 % verbose
            opts.qp.printLevel = 2; % iterations
        else
            opts.qp.printLevel = 0; % none
        end
    end
    %----------------------------------------------------------------------
    otherwise
        error('ERROR: unknown solver')
end
end