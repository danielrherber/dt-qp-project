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

% initialize quasilinearization-specific options structure
if ~isfield(opts,'qlin')
    opts.qlin = [];
end

% method to solve the nonlinear dynamic optimization problem
if ~isfield(opts.qlin,'method')
    opts.qlin.method = 'none';
    % opts.qlin.method = 'ipfmincon'; % fmincon nonlinear solver
    % opts.qlin.method = 'qlin'; % quasilinearization
end

% handle OLQ elements properly
if ~isfield(opts.qlin,'olqflag')
    % opts.qlin.olqflag = false; % everything nonlinear
    opts.qlin.olqflag = true; % direct incorporation of OLQ elements
end

%--------------------------------------------------------------------------
% ipfmincon options
%-------------------------------------------------------------------------
if strcmpi(opts.qlin.method,'ipfmincon')
    opts.qp.solver = 'ipfmincon';

    % optimization problem derivatives flag (NEED TO REMOVE)
    if ~isfield(opts.qlin,'derivativeflag')
        % opts.qlin.derivativeflag = false; % use finite differencing
        opts.qlin.derivativeflag = true; % use symbolic derivatives
    end

    % optimization problem derivatives method
    if ~isfield(opts.qlin,'derivativemethod')
        % opts.qlin.derivativemethod = 'internal'; % use internal fmincon real finite differencing
        % opts.qlin.derivativemethod = 'complex'; % use complex-step finite differencing
        opts.qlin.derivativemethod = 'symbolic'; % use symbolic derivatives
    end

end
%--------------------------------------------------------------------------
% quasilinearization options
%--------------------------------------------------------------------------
% relative function tolerance
if ~isfield(opts.qlin,'tolerance')
    opts.qlin.tolerance = 1e-6;
end

% maximum number of iterations
if ~isfield(opts.qlin,'imax')
    opts.qlin.imax = 50;
end

% constant scaling based on previous solution
if ~isfield(opts.qlin,'deltascaleflag')
    % opts.qlin.deltascaleflag = true; % enabled
    opts.qlin.deltascaleflag = false; % disabled
end

% sequential quadratic programming flag
if ~isfield(opts.qlin,'sqpflag')
    % opts.qlin.sqpflag = true; % enabled
    opts.qlin.sqpflag = false; % disabled
end
if opts.qlin.sqpflag
	opts.qlin.deltascaleflag = true; % required to be enabled

    % mirror negative eigenvalues in hessian
    if ~isfield(opts.qlin,'mirrorflag')
        % opts.qlin.mirrorflag = true; % enabled
        opts.qlin.mirrorflag = false; % disabled
    end

end

% trust region flag
if ~isfield(opts.qlin,'trustregionflag')
    % opts.qlin.trustregionflag = true; % enabled
    opts.qlin.trustregionflag = false; % disabled
end
if opts.qlin.trustregionflag
	opts.qlin.deltascaleflag = true; % required to be enabled

    % default trust region size
    if ~isfield(opts.qlin,'delta')
        opts.qlin.delta = 1;
    end

    % default contraction factor
    if ~isfield(opts.qlin,'xi')
        opts.qlin.xi = 0.8;
    end

end

% try to improve the initial guess value
if ~isfield(opts.qlin,'improveX0flag')
    % opts.qlin.improveX0flag = true; % enabled
    opts.qlin.improveX0flag = false; % disabled
end

end