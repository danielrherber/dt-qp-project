%--------------------------------------------------------------------------
% DTQP_createf.m
% Create the matrix for the linear objective function terms
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function f = DTQP_createf(l,m,p,opts)

    % initialize
	fJ = []; fV = [];

    % Lagrange terms
    if ~isempty(l)
        [~,J,V] = DTQP_L(l,p,opts);
        fJ = [fJ;J]; fV = [fV;V];
    end

    % Mayer terms
    if ~isempty(m)
        [~,J,V] = DTQP_M(m,p,opts);
        fJ = [fJ;J]; fV = [fV;V];
    end

    % sparse matrix for Hessian
    if isempty(fV)
        f = []; % no Hessian
    else
        % sparse matrix for gradient
        f = sparse(fJ,1,fV,p.nx,1);
    end

end