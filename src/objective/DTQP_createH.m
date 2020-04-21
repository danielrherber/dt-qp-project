%--------------------------------------------------------------------------
% DTQP_createH.m
% Create the matrix for the quadratic objective function terms
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function H = DTQP_createH(L,M,in,opts)

% initialize
HI = []; HJ = []; HV = [];

% Lagrange terms
if ~isempty(L)
    [I,J,V] = DTQP_L(L,in,opts);
    HI = [HI;I]; HJ = [HJ;J]; HV = [HV;V];
end

% Mayer terms
if ~isempty(M)
    [I,J,V] = DTQP_M(M,in,opts);
    HI = [HI;I]; HJ = [HJ;J]; HV = [HV;V];
end

% sparse matrix for Hessian
if isempty(HV)
    H = []; % no Hessian
else
    H = sparse(HI,HJ,HV,in.nx,in.nx);
    H = (H+H'); % make symmetric, then times 2 for 1/2*x'*H*x form
end

end