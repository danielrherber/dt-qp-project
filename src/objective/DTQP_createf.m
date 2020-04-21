%--------------------------------------------------------------------------
% DTQP_createf.m
% Create the matrix for the linear objective function terms
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function f = DTQP_createf(l,m,in,opts)

% initialize
fJ = []; fV = [];

% Lagrange terms
if ~isempty(l)
    [~,J,V] = DTQP_L(l,in,opts);
    fJ = [fJ;J]; fV = [fV;V];
end

% Mayer terms
if ~isempty(m)
    [~,J,V] = DTQP_M(m,in,opts);
    fJ = [fJ;J]; fV = [fV;V];
end

% sparse matrix for Hessian
if isempty(fV)
    f = []; % no Hessian
else
    % sparse matrix for gradient
    f = sparse(fJ,1,fV,in.nx,1);
end

end