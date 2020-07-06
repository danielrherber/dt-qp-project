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
    case 'qlin'
    % solve the NLDO problem with quasilinearization
    [T,U,Y,P,F,in,opts] = DTQP_QLIN(setup,opts);
    %------------------------------------------------------------------
    case 'nonlinearprogram'
    % solve the NLDO problem with interior point fmincon
    [T,U,Y,P,F,in,opts] = DTQP_IPFMINCON(setup,opts);
    %------------------------------------------------------------------
    otherwise
    error('opts.method.form not found')
    %------------------------------------------------------------------
end

end