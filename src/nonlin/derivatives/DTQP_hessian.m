%--------------------------------------------------------------------------
% DTQP_hessian.m
% Compute Hessian of the input function
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function D2ft = DTQP_hessian(in,auxdata,t,X,param,flag,k)

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
    D2ft = DTQP_tmultiprod(D2fi,auxdata,t);
    %----------------------------------------------------------------------
    case 'complex'
    D2ft = DTQP_hessian_complex_step(in.f{k},X,t,param);
    %----------------------------------------------------------------------
    case 'real-forward'
    D2ft = DTQP_hessian_real_forward(in.f{k},X,t,param);
    %----------------------------------------------------------------------
    case 'real-central'
    D2ft = DTQP_hessian_real_central(in.f{k},X,t,param);
    %----------------------------------------------------------------------
end

end