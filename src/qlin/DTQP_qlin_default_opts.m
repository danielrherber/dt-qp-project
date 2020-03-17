%--------------------------------------------------------------------------
% DTQP_qlin_default_opts.m
% Default options for the quasilinearization methods
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function opts = DTQP_qlin_default_opts(opts)

% relative function tolerance
if ~isfield(opts.qlin,'tolerance')
    opts.qlin.tolerance = 1e-6;
end

% maximum number of iterations
if ~isfield(opts.qlin,'imax')
    opts.qlin.imax = 500;
end

end