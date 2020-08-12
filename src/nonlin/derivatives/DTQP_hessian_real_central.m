%--------------------------------------------------------------------------
% DTQP_hessian_real_central.m
% Compute Hessian of the function using central real-step differentiation
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function D2f = DTQP_hessian_real_central(f,X,T,param)

% number of optimization variables
nx = size(X,2);

% number of time points
nt = length(T);

% initialize
D2f = zeros(nt,nx,nx);

% magnitude of the optimization variables
Xa = 1 + abs(X);

% differentiation step size
hi = (eps)^(1/3);
hj = (eps)^(1/3);

% step size vectors
Hi = hi*Xa;
Hj = hj*Xa;

% go through each optimization variable
for jx = 1:nx

    % loop for off diagonal Hessian
    for ix = jx:nx

        % compute multi-directional steps
        Xff = X;
        Xff(:,ix) = Xff(:,ix) + Hi(:,ix);
        Xff(:,jx) = Xff(:,jx) + Hj(:,jx);

        Xfb = X;
        Xfb(:,ix) = Xfb(:,ix) + Hi(:,ix);
        Xfb(:,jx) = Xfb(:,jx) - Hj(:,jx);

        Xbf = X;
        Xbf(:,ix) = Xbf(:,ix) - Hi(:,ix);
        Xbf(:,jx) = Xbf(:,jx) + Hj(:,jx);

        Xbb = X;
        Xbb(:,ix) = Xbb(:,ix) - Hi(:,ix);
        Xbb(:,jx) = Xbb(:,jx) - Hj(:,jx);

        % compute steps in both directions
        Uff = f(T,param,Xff);
        Ufb = f(T,param,Xfb);
        Ubf = f(T,param,Xbf);
        Ubb = f(T,param,Xbb);

        % compute second derivative
        D2f(:,jx,ix) = (Uff - Ufb - Ubf + Ubb)./(4*Hi(:,ix).*Hj(:,jx));

        % symmetric
        D2f(:,ix,jx) = D2f(:,jx,ix);

    end
end

end