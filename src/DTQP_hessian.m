%--------------------------------------------------------------------------
% DTQP_hessian.m
% Compute Hessian of the input function
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function D2ft = DTQP_hessian(in,p,t,X,param,flag,k)

% determine which derivative method
switch flag
    %----------------------------------------------------------------------
    case 'symbolic'
    % check if empty D2f (so all zeros)
    if isempty(in.D2f)
        D2ft = [];
        return
    end
    D2fi = DTQP_QLIN_update_tmatrix(in.D2f{k},[],X,param);
    D2ft = DTQP_tmultiprod(D2fi,p,t);
    %----------------------------------------------------------------------
    case 'complex'
    D2ft = DTQP_hessian_complex_step(in.f{k},X,t,param);
    %----------------------------------------------------------------------
    case 'real'
    error("Not implemented")
    %----------------------------------------------------------------------
end

end