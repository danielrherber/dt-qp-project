%--------------------------------------------------------------------------
% DTQP_createc.m
% Compute the constant objective term
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function c =  DTQP_createc(L,M,p,opts)

    % initialize constant term
	c = 0;            

    % define as a structure
    l.matrix = L;

    % add constant Lagrange term using same quadrature method as QP
    if ~isempty(L)
        [~,~,V] = DTQP_L(l,p,opts);
        c = sum(V); % calculate the integral
    end

    % add constant Mayer term
    if ~isempty(M)
        c = c + M;
    end

end