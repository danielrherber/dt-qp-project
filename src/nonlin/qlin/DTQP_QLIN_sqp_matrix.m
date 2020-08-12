%--------------------------------------------------------------------------
% DTQP_QLIN_sqp_matrix.m
% Construct sqp penalty matrix from sequences
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [H,Hsqp] = DTQP_QLIN_sqp_matrix(H,D2matrix,in,opts)

% initialize
HI = []; HJ = []; HV = [];

% create indices for Lagrangian penalty matrix
[I,J,V] = DTQP_SQP_lagrangian_penalty_matrix(D2matrix,in,opts);
HI = [HI;I]; HJ = [HJ;J]; HV = [HV;V];

% sparse matrix for Hessian
if isempty(HV)
    Hsqp = []; % no Hessian
else
    Hsqp = sparse(HI,HJ,HV,in.nx,in.nx);
    Hsqp = (Hsqp+Hsqp'); % make symmetric, then times 2 for 1/2*x'*H*x form
end

% combine with original hessian
if isempty(H)
    H = Hsqp;
else
    % combine
    H = H + Hsqp;

    % tolerance for symmetric positive semidefiniteness
    etol = sqrt(eps);

    if opts.method.mirrorflag

        % check if the matrix is symmetric positive semidefinite
        if eigs(H,1,'smallestreal') >= -etol
            % disp('Matrix is symmetric positive definite')
        else
            % disp('Matrix is not symmetric positive definite')

            % mirrored version
            [Um,Tm] = schur(full(H));
            H = Um*abs(Tm)*Um';
        end

        % make symmetric, then times 2 for 1/2*x'*H*x form
        H = (H+H')/2;
    end

end

end