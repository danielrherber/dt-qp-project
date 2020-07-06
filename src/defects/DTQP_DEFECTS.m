%--------------------------------------------------------------------------
% DTQP_DEFECTS.m
% Create the Aeq and beq matrices that represent the defect constraints
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [Aeq,beq,in] = DTQP_DEFECTS(A,B,G,d,in,opts)

switch upper(opts.dt.defects)
    case 'ZO' % zero-order hold
        [Aeq,beq] = DTQP_DEFECTS_ZO(A,B,G,d,in,opts);
    case 'EF' % Euler forward
        [Aeq,beq] = DTQP_DEFECTS_EF(A,B,G,d,in,opts);
    case 'TR' % trapezoidal
        [Aeq,beq] = DTQP_DEFECTS_TR(A,B,G,d,in,opts);
    case 'HS' % Hermite-Simpson
        [Aeq,beq] = DTQP_DEFECTS_HS(A,B,G,d,in,opts);
    case 'RK4' % fourth-order Runge-Kutta
        [Aeq,beq] = DTQP_DEFECTS_RK4(A,B,G,d,in,opts);
    case 'PS' % pseudospectral (both LGL and CGL)
        [Aeq,beq] = DTQP_DEFECTS_PS(A,B,G,d,in,opts);
    case 'HUEN' % Heun's method
        [Aeq,beq] = DTQP_DEFECTS_Huen(A,B,G,d,in,opts);
    case 'MODEF' % Modified Euler method
        [Aeq,beq] = DTQP_DEFECTS_ModEF(A,B,G,d,in,opts);
end

% (potentially) store defect constraint indices
if isfield(opts.method,'sqpflag') && opts.method.sqpflag
    in.multipliers.defects = reshape(1:size(Aeq,1),[],in.ny);
end

end