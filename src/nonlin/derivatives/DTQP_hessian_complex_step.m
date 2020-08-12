%--------------------------------------------------------------------------
% DTQP_hessian_complex_step.m
% Compute Hessian of the function using complex-step differentiation
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function D2f = DTQP_hessian_complex_step(f,X,T,param)

% number of optimization variables
nx = size(X,2);

% number of time points
nt = length(T);

% initialize
D2f = zeros(nt,nx,nx);

% differentiation step size
hj = eps;
hi = (eps)^(1/3);

% square of step size
h2 = hj*hi*2;

% go through each optimization variable
for jx = 1:nx

    % reference point
    X0 = X;

    % increment in independent variable
    X0(:,jx) = X0(:,jx) + hj*1i;

    % loop for off diagonal Hessian
    for ix = jx:nx

        % reference with increment
        Xp = X0;
        Xm = X0;

        % positive real increment
        Xp(:,ix) = Xp(:,ix) + hi;

        % function call with a double increment
        up = f(T,param,Xp);

        % negative real increment
        Xm(:,ix) = Xm(:,ix) - hi;

        % function call with a double increment
        um = f(T,param,Xm);

        % Hessian (central + complex step)
        D2f(:,jx,ix) = imag(up-um)/h2;

        % symmetric
        D2f(:,ix,jx) = D2f(:,jx,ix);

    end

end