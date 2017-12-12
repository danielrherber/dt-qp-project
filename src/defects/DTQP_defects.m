%--------------------------------------------------------------------------
% DTQP_defects.m
% Create the Aeq and beq matrices that represent the defect constraints
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [Aeq,beq] = DTQP_defects(A,B,G,d,p,opts)

    switch opts.Defectmethod
        case 'ZO' % zero-order hold
            [Aeq,beq] = DTQP_defects_ZO(A,B,G,d,p,opts);
        case 'EF' % Euler forward 
            [Aeq,beq] = DTQP_defects_EF(A,B,G,d,p,opts);
        case 'TR' % trapezoidal
            [Aeq,beq] = DTQP_defects_TR(A,B,G,d,p,opts);
        case 'HS' % Hermite-Simpson 
            [Aeq,beq] = DTQP_defects_HS(A,B,G,d,p,opts); 
        case 'RK4' % fourth-order Runge-Kutta 
            [Aeq,beq] = DTQP_defects_RK4(A,B,G,d,p,opts);   
        case 'PS' % pseudospectral (both LGL and CGL)
            [Aeq,beq] = DTQP_defects_PS(A,B,G,d,p,opts);
    end
    
end