%--------------------------------------------------------------------------
% DTQP_solve.m
% Construct and solve the LQDO problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts)

    % initialize some stuff
    [setup,opts] = DTQP_default_opts(setup,opts);
    
    % potentially start the timer
    if (opts.general.displevel > 0) % minimal
        tic % start timer
    end
    
    % check for a multiphase problem
    if length(setup) > 1
        solvefun = @DTQP_multiphase;
    else
        solvefun = @DTQP_singlephase;
    end

    % solve the problem potentially using mesh refinement 
    [T,U,Y,P,F,in,opts] = DTQP_meshr(setup,opts,solvefun);

end