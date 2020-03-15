%--------------------------------------------------------------------------
% DTQP_meshr.m
% Solve the LQDO problem with the selected mesh refinement scheme
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [T,U,Y,P,F,in,opts] = DTQP_meshr(setup,opts)

% determine which mesh refinement scheme should be used
switch upper(opts.dt(1).meshr.method)
    %----------------------------------------------------------------------
    case 'NONE' % no mesh refinement
        [T,U,Y,P,F,in,opts] = DTQP_multiphase(setup,opts);

        % end the timer
        if (opts.general.displevel > 0) % minimal
            in(end).QPtotaltime = toc;
        end

        % display to the command window
        if (opts.general.displevel > 1) % verbose
            disp('----------------------------------------')
            disp(['QP total time: ', num2str(in(end).QPtotaltime), ' s'])
        end
    %----------------------------------------------------------------------
    case 'RICHARDSON-DOUBLING' % doubling the mesh and Richardson extrapolation
        [T,U,Y,P,F,in,opts] = DTQP_meshr_richardson_doubling(setup,opts);
    %----------------------------------------------------------------------
    % case 'CUBIC-BSPLINES' % in development
    %     [T,U,Y,P,F,in,opts] = DTQP_meshr_cubic_bsplines(setup,opts,solvefun);
    %----------------------------------------------------------------------
    case 'TEST' % in development
        [T,U,Y,P,F,in,opts] = DTQP_meshr_test(setup,opts);
    %----------------------------------------------------------------------
    otherwise
        error('mesh refinement method invalid')
    %----------------------------------------------------------------------
end

end