%--------------------------------------------------------------------------
% DTQP_NLP_parse_additional_constraints.m
% Parse element.h or element.g to obtain nonlinear constraint function and
% linear constraint matrices (as needed)
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Project Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [A,b,Ic,in,opts,setup] = DTQP_NLP_parse_additional_constraints(setup,in,opts,olqflag,nt,nI,cflag)

% determine which type of constraint
if cflag % inequality constraints

    % check if there is a g field
    if isfield(setup.element,'g')

        % extract
        c = setup.element.g;
        YZ_original = setup.Z;
        cstr = 'cin';

    else

        % create linear path and boundary inequality constraints
        [A,b] = DTQP_create_YZ(setup.Z,in);

        % default values
        in.cin = [];
        Ic = [];

        return
    end

else % equality constraints

    % check if there is a h field
    if isfield(setup.element,'h')

        % extract
        c = setup.element.h;
        YZ_original = setup.Y;
        cstr = 'ceq';

    else % only linear equality constraints

        % create linear path and boundary equality constraints
        [A,b] = DTQP_create_YZ(setup.Y,in);

        % default values
        in.ceq = [];
        Ic = [];

        return

    end

end

% check if pathboundary was provided
if isfield(c,'pathboundary')
    pathboundary = c.pathboundary;
end

% check internal information is already available for the additional constraints
if isfield(in.internalinfo,cstr)

    % extract
    if cflag
        c = in.internalinfo.cin;
    else
        c = in.internalinfo.ceq;
    end

else

    % calculate derivatives
    [c,opts] = DTQP_NLP_symb(c.func,in,olqflag,true,opts);

    % store to internal information structure
    if cflag
        in.internalinfo.cin = c;
    else
        in.internalinfo.ceq = c;
    end

end

% number of constraints
nz = length(c.f);

% assign pathboundary or extract if determined symbolically
if isfield(c,'pathboundary')
	pathboundary = c.pathboundary;
else
    c.pathboundary = pathboundary;
end

% initialize
Ic = cell(nz,1);

% go through each constraint
for k = 1:nz

    % check if this is a path constraint
    if pathboundary(k)
        Ii = nI + (1:nt)';
    else
        Ii = nI + 1;
    end

    % update
    nI = max(Ii);

    % assign
    Ic{k} = Ii;

end

% store to internal information structure
if cflag
    in.internalinfo.cin = c;
    in.cin = c;
    in.Ilambda.cin = Ic;
else
    in.internalinfo.ceq = c;
    in.ceq = c;
    in.Ilambda.ceq = Ic;
end

% check if there are any linear constraints to add
nnew = length(c.Ilin);

% initialize structure
if nnew > 0
    YZ = struct('linear',[],'b',[]);
    YZ(nnew) = YZ;
else
    YZ = [];
end

% create additional YZ structure
for k = 1:nnew

    % reset structure and local index
    YZ_ = [];
    idx = 0;

    % controls
    if ~isempty(c.B)

        % increment local index
        idx = idx + 1;

        % variable
        YZ_.linear(idx).right = 1; % controls

        % assign
        YZ_.linear(idx).matrix = c.B(k,:);

    end

    % states
    if ~isempty(c.A)

        % increment local index
        idx = idx + 1;

        % variable
        YZ_.linear(idx).right = 2; % states

        % assign
        YZ_.linear(idx).matrix = c.A(k,:);

    end

    % parameters
    if ~isempty(c.G)

        % increment local index
        idx = idx + 1;

        % variable
        YZ_.linear(idx).right = 3; % parameters

        % assign
        YZ_.linear(idx).matrix = c.G(k,:);

    end

    % initial states
    if ~isempty(c.Ai)

        % increment local index
        idx = idx + 1;

        % variable
        YZ_.linear(idx).right = 4; % initial states

        % assign
        YZ_.linear(idx).matrix = c.Ai(k,:);

    end

    % final states
    if ~isempty(c.Af)

        % increment local index
        idx = idx + 1;

        % variable
        YZ_.linear(idx).right = 5; % final states

        % assign
        YZ_.linear(idx).matrix = c.Af(k,:);

    end

    % constant (negative disturbance)
    if ~isempty(c.d)

        % assign
        YZ_.b = {'prod',-1,c.d(k,:)}; % constant

    end

    % assign
    YZ(k) = YZ_;

end

% combine with previously defined linear constraints
YZ = [YZ_original;YZ];

% create linear path and boundary constraints
[A,b] = DTQP_create_YZ(YZ,in);

% assign (note: might not be needed anymore)
if cflag
    setup.Z = YZ;
else
    setup.Y = YZ;
end

end