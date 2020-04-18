%--------------------------------------------------------------------------
% DTQP_qlin_sqpMatrix.m
% Construct sqp penalty matrix from sequences
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function H = DTQP_qlin_sqpMatrix(D2matrix,in,opts)

% initialize
HI = []; HJ = []; HV = [];

% create indices for Lagrangian penalty matrix
[I,J,V] = DTQP_SQP_lagrangianPenaltyMatrix(D2matrix,in,opts);
HI = [HI;I]; HJ = [HJ;J]; HV = [HV;V];

% sparse matrix for Hessian
if isempty(HV)
    H = []; % no Hessian
else
    H = sparse(HI,HJ,HV,in.nx,in.nx);
    H = (H+H'); % make symmetric, then times 2 for 1/2*x'*H*x form
end

end