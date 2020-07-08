%--------------------------------------------------------------------------
% DTQP_MESH.m
% Solve the LQDO problem with the selected mesh refinement scheme
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [T,U,Y,P,F,in,opts] = DTQP_MESH(setup,opts)

% determine which mesh refinement scheme should be used
switch upper(opts.dt(1).meshr.method)
    %----------------------------------------------------------------------
    case 'NONE' % no mesh refinement
        [T,U,Y,P,F,in,opts] = DTQP_multiphase(setup,opts);
    %----------------------------------------------------------------------
    case 'RICHARDSON-DOUBLING' % doubling the mesh and Richardson extrapolation
        [T,U,Y,P,F,in,opts] = DTQP_MESH_richardson_doubling(setup,opts);
    %----------------------------------------------------------------------
    case 'SS-BETTS' % single-step method mesh refinement from Betts textbook
        [T,U,Y,P,F,in,opts] = DTQP_MESH_ss_betts(setup,opts);
    %----------------------------------------------------------------------
    case 'TEST' % in development
        % [T,U,Y,P,F,in,opts] = DTQP_meshr_test(setup,opts);
    %----------------------------------------------------------------------
    otherwise
        error('mesh refinement method invalid')
    %----------------------------------------------------------------------
end

end