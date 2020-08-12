%--------------------------------------------------------------------------
% DTQP_jacobian.m
% Compute Jacobian of the input functions
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function Dft = DTQP_jacobian(in,p,t,X,param,flag)

% determine which derivative method
switch flag
    %----------------------------------------------------------------------
    case 'symbolic'
    Dfi = DTQP_QLIN_update_tmatrix(in.Df,[],X,param);
    Dft = DTQP_tmultiprod(Dfi,p,t);
    %----------------------------------------------------------------------
    case 'complex'
    Dft = DTQP_jacobian_complex_step(in.f,X,t,param);
    %----------------------------------------------------------------------
    case 'real-forward'
    Dft = DTQP_jacobian_real_forward(in.f,X,t,param);
    %----------------------------------------------------------------------
    case 'real-central'
    Dft = DTQP_jacobian_real_central(in.f,X,t,param);
    %----------------------------------------------------------------------
end

end