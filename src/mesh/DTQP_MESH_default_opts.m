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
    % meshr.method = 'cubic-bsplines'; % in development
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
    % case 'CUBIC-BSPLINES' % see DTQP_meshr_cubic_bsplines.m
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
end

end