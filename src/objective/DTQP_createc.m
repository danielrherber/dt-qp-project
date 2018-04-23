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
function c =  DTQP_createc(cL,cM,p,opts)

    % initialize constant term
	c = 0;            

    % add constant Lagrange term using same quadrature method as QP
    if ~isempty(cL)
        [~,~,V] = DTQP_L(cL,p,opts);
        c = sum(V); % calculate the integral
    end

    % add constant Mayer term
    if ~isempty(cM)
        c = c + cM;
    end

end