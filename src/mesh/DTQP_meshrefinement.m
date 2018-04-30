%--------------------------------------------------------------------------
% DTQP_meshrefinement.m
% Solve the LQDO problem with the selected mesh refinement strategy
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [T,U,Y,P,F,p,opts] = DTQP_meshrefinement(setup,opts,solvefun)

% determine which mesh refinement scheme should be used
switch lower(opts.dt(1).refinement)
    %----------------------------------------------------------------------
    case 'none' % no mesh refinement
        [T,U,Y,P,F,p,opts] = solvefun(setup,opts);
    %---------------------------------------------------------------------- 
    otherwise
        error('mesh refinement method invalid')
    %----------------------------------------------------------------------
end

end