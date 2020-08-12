%--------------------------------------------------------------------------
% DTQP_hessian_real_forward.m
% Compute Hessian of the function using forward real-step differentiation
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function D2f = DTQP_hessian_real_forward(f,X,T,param)

% number of optimization variables
nx = size(X,2);

% number of time points
nt = length(T);

% initialize
D2f = zeros(nt,nx,nx);

% differentiation step size
hj = (eps)^(1/3);
hi = (eps)^(1/3);

% product of step sizes
hij = hi*hj;

% reference point matrix
X0 = repmat(X,[nx 1]);

% function call at original point
U0 = f(T,param,X);

% replicate parameter vector if it is time-dependent
if size(param,1) == nt
    PARAM = repmat(param,[nx 1]);
else
    PARAM = param;
end

% create expanded step size matrices in first direction
Hi = zeros(nx*nt,nx);
for ix = 1:nx
    Hi(((ix-1)*nt+1):(ix*nt),ix) = hi;
end

% compute forward steps in first direction
Ui = f(T,PARAM,X0 + Hi);
Ui = reshape(Ui,[],nx);

% check if step sizes are the same
if hi == hj

    % use previous data
    Uj = Ui;

else

    % create expanded step size matrices in second direction
    Hj = zeros(nx*nt,nx);
    for jx = 1:nx
        Hj(((jx-1)*nt+1):(jx*nt),jx) = hj;
    end

    % compute forward steps in second direction
    Uj = f(T,PARAM,X0 + Hj);
    Uj = reshape(Uj,[],nx);

end

% unit step size vectors
ei = diag(repelem(hi,nx));
ej = diag(repelem(hj,nx));

% go through each optimization variable
for jx = 1:nx

    % loop for off diagonal Hessian
    for ix = jx:nx

        % compute multi-directional step
        Xij = X + ei(ix,:) + ej(jx,:);

        % compute forward steps in both directions
        Uij = f(T,param,Xij);

        % compute second derivative
        D2f(:,jx,ix) = (Uij - Ui(:,ix) - Uj(:,jx) + U0)/hij;

        % symmetric
        D2f(:,ix,jx) = D2f(:,jx,ix);

    end
end

end