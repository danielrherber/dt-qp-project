%--------------------------------------------------------------------------
% DTQP_solve.m
% Construct and solve the LQDO problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [T,U,Y,P,F,p,opts] = DTQP_solve(setup,opts)

    % initialize some stuff
    [setup,opts] = DTQP_default_opts(setup,opts);

    % check for a multiphase problem
    if length(setup) > 1
        solvefun = @DTQP_multiphase;
    else
        solvefun = @DTQP_singlephase;
    end

    % solve the problem potentially using mesh refinement 
    [T,U,Y,P,F,p,opts] = DTQP_meshrefinement(setup,opts,solvefun);

end