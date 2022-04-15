%--------------------------------------------------------------------------
% DTQP_NONLIN.m
% Solve the nonlinear DO problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [T,U,Y,P,F,in,opts] = DTQP_NONLIN(setup,opts)

% determine which nonlinear method should be used
switch opts.method.form
    %------------------------------------------------------------------
    case {'nonlinearprogram','nlp'}
    % solve the NLDO problem with nonlinear programming
    [T,U,Y,P,F,in,opts] = DTQP_NLP(setup,opts);
    %------------------------------------------------------------------
    case 'qlin'
    % solve the NLDO problem with quasilinearization
    [T,U,Y,P,F,in,opts] = DTQP_QLIN(setup,opts);
    %------------------------------------------------------------------
    otherwise
    error('opts.method.form not found')
    %------------------------------------------------------------------
end

end