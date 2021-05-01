%--------------------------------------------------------------------------
% DTQP_MESH_default_opts.m
% Default options for the mesh refinement schemes
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function meshr = DTQP_MESH_default_opts(meshr)

% mesh refinement method
if ~isfield(meshr,'method')
    meshr.method = 'none'; % no mesh refinement
    % meshr.method = 'richardson-doubling'; % Richardson doubling method
    % meshr.method = 'ss-betts';
end

switch upper(meshr.method)
    %----------------------------------------------------------------------
    case 'NONE' % see DTQP_solve.m

    % no options needed

    %----------------------------------------------------------------------
    case 'RICHARDSON-DOUBLING' % see DTQP_MESH_richardson_doubling.m

    % relative function tolerance
    if ~isfield(meshr,'tolerance')
        meshr.tolerance = 1e-5;
    end

    % maximum number of time points
    if ~isfield(meshr,'ntmax')
        meshr.ntmax = 30000;
    end

    % initial number of time points (only used if nt = nan)
    if ~isfield(meshr,'ntinit')
        meshr.ntinit = 5;
    end

    %----------------------------------------------------------------------
    case 'SS-BETTS' % see DTQP_MESH_ss_betts.m

    % tolerance for the maximum relative local error
    if ~isfield(meshr,'tolerance')
        meshr.tolerance = 1e-4;
    end

    % predicted error safety factor (0 < errorsafety < 1)
    if ~isfield(meshr,'kappa')
        meshr.errorsafety = 1/10;
    end

    % maximum number of points to add to a single interval
    if ~isfield(meshr,'ntmaxinterval')
        meshr.ntmaxinterval = 5;
    end

    % number of integration points when computing absolute local error
    if ~isfield(meshr,'ntintegration')
        meshr.ntintegration = 10;
    end

    % maximum number of mesh iterations
    if ~isfield(meshr,'maxiters')
        meshr.maxiters = 20;
    end

    % maximum number of time points
    if ~isfield(meshr,'ntmax')
        meshr.ntmax = 10000;
    end

    % flag to save the intermediate meshes (with errors)
    if ~isfield(meshr,'savemesh')
        % meshr.storemesh = false; % don't save the intermediate meshes
        meshr.storemesh = true; % save the intermediate meshes
    end

    %----------------------------------------------------------------------
    case 'TEST' % see .m

    % relative function tolerance
    if ~isfield(meshr,'tolerance')
        meshr.tolerance = 1e-5;
    end

    % maximum number of time points
    if ~isfield(meshr,'ntmax')
        meshr.ntmax = 30000;
    end

    % initial number of time points (only used if nt = nan)
    if ~isfield(meshr,'ntinit')
        meshr.ntinit = 5;
    end

    %----------------------------------------------------------------------
    otherwise

    error('ERROR: unknown mesh refinement scheme')

    %----------------------------------------------------------------------
end

end