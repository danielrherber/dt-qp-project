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
    % check if there are symbolically defined problem elements
    if isfield(setup,'element')
        solve_fun = @DTQP_NONLIN; % NLDO problem
    else
        solve_fun = @DTQP_multiphase; % LQDO problem
    end

    % solve the problem
    [T,U,Y,P,F,in_,opts] = solve_fun(setup,opts);
    %----------------------------------------------------------------------
    case 'RICHARDSON-DOUBLING' % doubling the mesh and Richardson extrapolation
    [T,U,Y,P,F,in_,opts] = DTQP_MESH_richardson_doubling(setup,opts);
    %----------------------------------------------------------------------
    case 'SS-BETTS' % single-step method mesh refinement from Betts textbook
    [T,U,Y,P,F,in_,opts] = DTQP_MESH_ss_betts(setup,opts);
    %----------------------------------------------------------------------
    case 'TEST' % in development
    % [T,U,Y,P,F,in,opts] = DTQP_meshr_test(setup,opts);
    %----------------------------------------------------------------------
    otherwise
    error('mesh refinement method invalid')
    %----------------------------------------------------------------------
end

% modify in structure
in.phase_info = in_;
in.output = in_(1).output;
in.auxdata = in_(1).auxdata;
in.t0 = in_(1).t0;
in.tf = in_(end).tf;
in.nt = sum([in_.nt]);

end