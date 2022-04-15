%--------------------------------------------------------------------------
% DTQP_NONLIN_default_opts.m
% Default options for the NLDO methods
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function opts = DTQP_NONLIN_default_opts(opts)

% method to solve the nonlinear dynamic optimization problem
if ~isfield(opts.method,'form')
    opts.method.form = 'none';
    % opts.method.form = 'nonlinearprogram'; % nonlinear programming
    % opts.method.form = 'qlin'; % quasilinearization
end

% handle OLQ elements properly
if ~isfield(opts.method,'olqflag')
    % opts.method.olqflag = false; % everything nonlinear
    opts.method.olqflag = true; % direct incorporation of OLQ elements
end

%--------------------------------------------------------------------------
% ipfmincon options
%-------------------------------------------------------------------------
if any(strcmpi(opts.method.form,{'nonlinearprogram','nlp'}))

    % set solver function
    opts.solver.function = 'ipfmincon';

    % optimization problem derivatives method
    if ~isfield(opts.method,'derivatives')
        % opts.method.derivatives = 'internal'; % use internal fmincon real finite differencing
        % opts.method.derivatives = 'real-forward'; % use forward real-step finite differencing
        % opts.method.derivatives = 'real-central'; % use central real-step finite differencing
        % opts.method.derivatives = 'complex'; % use forward complex-step finite differencing
        opts.method.derivatives = 'symbolic'; % use symbolic derivatives
    end

end
%--------------------------------------------------------------------------
% quasilinearization options
%--------------------------------------------------------------------------
if strcmpi(opts.method.form,'qlin')

    % relative function tolerance
    if ~isfield(opts.method,'tolerance')
        opts.method.tolerance = 1e-6;
    end

    % maximum number of iterations
    if ~isfield(opts.method,'maxiters')
        opts.method.maxiters = 50;
    end

    % constant scaling based on previous solution
    if ~isfield(opts.method,'deltascaleflag')
        % opts.method.deltascaleflag = true; % enabled
        opts.method.deltascaleflag = false; % disabled
    end

    % sequential quadratic programming flag
    if ~isfield(opts.method,'sqpflag')
        % opts.method.sqpflag = true; % enabled
        opts.method.sqpflag = false; % disabled
    end
    if opts.method.sqpflag
        opts.method.deltascaleflag = true; % required to be enabled

        % mirror negative eigenvalues in hessian
        if ~isfield(opts.method,'mirrorflag')
            % opts.method.mirrorflag = true; % enabled
            opts.method.mirrorflag = false; % disabled
        end

    end

    % trust region flag
    if ~isfield(opts.method,'trustregionflag')
        % opts.method.trustregionflag = true; % enabled
        opts.method.trustregionflag = false; % disabled
    end
    if opts.method.trustregionflag
        opts.method.deltascaleflag = true; % required to be enabled

        % default trust region size
        if ~isfield(opts.method,'delta')
            opts.method.delta = 1;
        end

        % default contraction factor
        if ~isfield(opts.method,'xi')
            opts.method.xi = 0.8;
        end

    end

    % try to improve the initial guess value
    if ~isfield(opts.method,'improveguess')
        % opts.method.improveguess = true; % enabled
        opts.method.improveguess = false; % disabled
    end

end

end