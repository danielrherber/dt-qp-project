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
switch opts.qlin.method
    %------------------------------------------------------------------
    case 'qlin'
    % solve the NLDO problem with quasilinearization
    [T,U,Y,P,F,in,opts] = DTQP_QLIN(setup,opts);
    %------------------------------------------------------------------
    case 'ipfmincon'
    % solve the NLDO problem with interior point fmincon
    [T,U,Y,P,F,in,opts] = DTQP_IPFMINCON(setup,opts);
    %------------------------------------------------------------------
    otherwise
    error('opts.qlin.method not found')
    %------------------------------------------------------------------
end

end