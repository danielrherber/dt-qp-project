%--------------------------------------------------------------------------
% DTQP_SOLVER_default_opts.m
% Defaults options for the optimization solvers
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function opts = DTQP_SOLVER_default_opts(opts)

switch upper(opts.solver.function)
    %----------------------------------------------------------------------
    case {'QUADPROG','BUILT-IN'} % see DTQP_SOLVER_quadprog.m
    % tolerance
    if ~isfield(opts.solver,'tolerance')
        opts.solver.tolerance = 1e-12;
    end

    % maximum iterations
    if ~isfield(opts.solver,'maxiters')
        opts.solver.maxiters = 200;
    end

    % display level in the optimization routine
    if ~isfield(opts.solver,'display')
        opts.solver.display = 'none'; % none
        if opts.general.displevel > 1 % verbose
            opts.solver.display = 'iter'; % iterations
        end
        if ~strcmpi(opts.dt(1).meshr.method,'none')
            opts.solver.display = 'none'; % iterations
        end
    end
    %----------------------------------------------------------------------
    case 'CVX' % see DTQP_SOLVER_cvx.m
    %----------------------------------------------------------------------
    case 'QPOASES' % see DTQP_SOLVER_qpoases.m
    % maximum iterations
    if ~isfield(opts.solver,'maxiters')
        opts.solver.maxiters = -1; % value chosen heuristically
    end

    % maximum CPU time in seconds
    if ~isfield(opts.solver,'maxcputime')
        opts.solver.maxcputime = -1; % only iteration limit is used
    end

    % display level in the optimization routine
    if ~isfield(opts.solver,'printLevel')
        if opts.general.displevel > 1 % verbose
            opts.solver.printLevel = 2; % iterations
        else
            opts.solver.printLevel = 0; % none
        end
    end
    %----------------------------------------------------------------------
    case 'OSQP' % see DTQP_SOLVER_osqp.m
    % tolerance
    if ~isfield(opts.solver,'tolerance')
        opts.solver.tolerance = 1e-4;
    end

    % maximum iterations
    if ~isfield(opts.solver,'maxiters')
        opts.solver.maxiters = 100000;
    end
    %----------------------------------------------------------------------
    case 'IPFMINCON' % see DTQP_SOLVER_ipfmincon.m
    % tolerance
    if ~isfield(opts.solver,'tolerance')
        opts.solver.tolerance = 1e-8;
    end

    % maximum iterations
    if ~isfield(opts.solver,'maxiters')
        opts.solver.maxiters = 200;
    end

    % display level in the optimization routine
    if ~isfield(opts.solver,'display')
        opts.solver.display = 'none'; % none
        if opts.general.displevel > 1 % verbose
            opts.solver.display = 'iter'; % iterations
        end
        if ~strcmpi(opts.dt(1).meshr.method,'none')
            opts.solver.display = 'none'; % iterations
        end
    end
    %----------------------------------------------------------------------
    otherwise
    error('ERROR: unknown solver')
    %----------------------------------------------------------------------
end
end