%--------------------------------------------------------------------------
% DTQP_NLP_constraints.m
% Compute nonlinear constraint values and Jacobians
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [G,H,DG,DH] = DTQP_NLP_constraints(X,dyn,cin,ceq,in,opts)

% check if gradients are requested
if nargout > 2
    Dflag = true;
else
    Dflag = false;
end

% compute defect constraints
% TODO: add more methods and move outside
if ~isfield(dyn,'Inon')
    dyn.Inon = 1:length(dyn.f);
end

% only if nonlinear state derivative functions
if ~isempty(dyn.Inon)
    switch upper(opts.dt.defects)
        case 'ZO' % zero-order hold
            error(' ')
        case 'EF' % Euler forward
            error(' ')
        case 'TR' % trapezoidal
            [z,Dz] = DTQP_DEFECTS_TR_nonlin(X,dyn,in,opts,Dflag);
        case 'HS' % Hermite-Simpson
%             error(' ')
            [z,Dz] = DTQP_DEFECTS_HS_nonlin(X,dyn,in,opts,Dflag);
        case 'RK4' % fourth-order Runge-Kutta
            error(' ')
        case 'PS' % pseudospectral (both LGL and CGL)
            [z,Dz] = DTQP_DEFECTS_PS_nonlin(X,dyn,in,opts,Dflag);
        case 'HUEN' % Heun's method
            error(' ')
        case 'MODEF' % Modified Euler method
            error(' ')
    end
else
    z = []; Dz = [];
end

% compute additional nonlinear equality constraints
if ~isempty(ceq)
    [h,Dh] = DTQP_NLP_additional_constraints(X,ceq,in,opts,Dflag);
else
    h = []; Dh = [];
end

% combine equality constraints
H = [z;h];

% compute additional nonlinear inequality constraints
if ~isempty(cin)
    [G,DG] = DTQP_NLP_additional_constraints(X,cin,in,opts,Dflag);
else
    G = []; DG = [];
end

% compute Jacobians if needed
if Dflag

    % output Jacobian of nonlinear equality constraints
    DH = [Dz;Dh]';

    % output Jacobian of nonlinear inequality constraints
    DG = DG';

end

end