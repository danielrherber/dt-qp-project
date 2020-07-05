%--------------------------------------------------------------------------
% DTQP_jacobian_complex_step.m
% Compute Jacobian of the functions using complex-step differentiation
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function Df = DTQP_jacobian_complex_step(f,X,T,param)

% number of functions
nf = length(f);

% number of optimization variables
nx = size(X,2);

% number of time points
nt = length(T);

% initialize
Df = zeros(nt,nf,nx);

% differentiation step size
h = nx*eps;

% go through each function
for kx = 1:nf

    % extract current function
    fun = f{kx};

    % go through each optimization variable
    for jx = 1:nx

        % reference point
        X0 = X;

        % increment optimization variable at each time point
        X0(:,jx) = X0(:,jx) + h*1i;

        % complex step differentiation
        Df(:,kx,jx) = imag(fun(T,param,X0))/h;

    end
end

end