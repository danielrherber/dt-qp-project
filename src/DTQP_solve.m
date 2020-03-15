%--------------------------------------------------------------------------
% DTQP_solve.m
% Construct and solve a (LQDO) problem with DTQP
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

% solve the problem potentially using mesh refinement
[T,U,Y,P,F,in,opts] = DTQP_meshr(setup,opts);

end