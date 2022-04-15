%--------------------------------------------------------------------------
% DTQP_NLP_parse_dynamics.m
% Parse element.dynamics to obtain nonlinear state derivative function and
% linear constraint matrices (as needed)
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Project Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [Aeq1,beq1,Idyn,in,opts,nI] = DTQP_NLP_parse_dynamics(setup,in,opts,olqflag,nt)

% check if there is a dynamics field
if isfield(setup.element,'dynamics')

    % extract
    D = setup.element.dynamics;

else % only linear dynamic equations

    % construct linear defect constraint matrices
    [Aeq1,beq1,in] = DTQP_DEFECTS(setup.A,setup.B,setup.G,setup.d,in,opts);

    % default values
    in.dyn = [];
    Idyn = [];
    nI = 0;

    return

end

% check internal information is already available for the dynamics
if isfield(in.internalinfo,'dyn')

    % extract
    dyn = in.internalinfo.dyn;

else

    % calculate derivatives
    [dyn,opts] = DTQP_NLP_symb(D,in,olqflag,false,opts);

    % store to internal information structure
    in.internalinfo.dyn = dyn;

end

% number of constraints
nz = length(dyn.f);

% initialize
Idyn = cell(nz,1);
nI = 0;

% go through each constraint
for k = 1:nz

    % check the method
    if strcmpi(opts.dt.defects,'PS')
        Ii = nI + (1:nt)';
    else
        Ii = nI + (1:nt-1)';
    end

    % update
    nI = max(Ii);

    % assign
    Idyn{k} = Ii;

end

% assign
in.dyn = dyn;
in.Ilambda.dyn = Idyn;

% current number of nonlinear equality constraints
nI = max(vertcat(in.Ilambda.dyn{:}));

% initialize linear defect constraint matrices
Aeq1 = zeros(0,in.nx); beq1 = zeros(0,1);

% check if this field is present
if isfield(dyn,'Ilin')

    % check if there are any linear dynamic equations
    if ~isempty(dyn.Ilin)

        % convert to DTQP compatible functions
        if ~isempty(dyn.A)
            setup.A = DTQP_QLIN_update_tmatrix(dyn.A,[],[],in.param);
        else
            setup.A = [];
        end
        if ~isempty(dyn.B)
            setup.B = DTQP_QLIN_update_tmatrix(dyn.B,[],[],in.param);
        else
            setup.B = [];
        end
        if ~isempty(dyn.G)
            setup.G = DTQP_QLIN_update_tmatrix(dyn.G,[],[],in.param);
        else
            setup.G = [];
        end
        if ~isempty(dyn.d)
            setup.d = DTQP_QLIN_update_tmatrix(dyn.d,[],[],in.param);
            in.nd = size(setup.d,1); % update number of disturbances
        else
            setup.d = [];
        end

        % assign
        in.IDlin = dyn.Ilin;
        in.IDnon = dyn.Inon;

        % construct linear defect constraints
        [Aeq1,beq1,in] = DTQP_DEFECTS(setup.A,setup.B,setup.G,setup.d,in,opts);

    end
end

end