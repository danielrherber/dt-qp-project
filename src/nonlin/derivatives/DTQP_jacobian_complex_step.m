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
h = eps;

% create expanded step size matrix
H2 = zeros(nx*nt,nx);
for jx = 1:nx
    H2(((jx-1)*nt+1):(jx*nt),jx) = h;
end

% reference point matrix
X0 = repmat(X,[nx 1]);

% compute forward step matrix
Xf = X0 + H2*1i;

% replicate parameter vector if it is time-dependent
if size(param,1) == nt
    PARAM = repmat(param,[nx 1]);
else
    PARAM = param;
end

% go through each function
for kx = 1:nf

    % extract current function
    fun = f{kx};

    % function call with forward increment
    Uf = fun(T,PARAM,Xf);

    % forward-step differentiation
    Dft = imag(Uf)/h;

    % reshape and place in Df
    Df(:,kx,:) = reshape(Dft,nt,1,nx);

end

end