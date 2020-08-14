%--------------------------------------------------------------------------
% DTQP_IPFMINCON.m
% Prepare and solve NLDO problem using interior point fmincon
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [T,U,Y,P,F,in,opts] = DTQP_IPFMINCON(setup,opts)

% initialize some stuff
[setup,in] = DTQP_initialize(setup,opts.dt);

% extract
nu = in.nu; ny = in.ny; np = in.np; nt = in.nt;
symb = setup.symb;

% initialize
lqdoflag = true;
in.Ilambda = [];
linflag = opts.method.olqflag;

% set default field value
if isfield(symb,'param')

    % convert time-varying cell to matrix on specified mesh
    % NOTE: this won't work for methods that need the parameters and
    % different points than the original mesh
    if isa(symb.param,'cell')
        symb.param = squeeze(DTQP_tmatrix(symb.param,setup.p,in.t));
    end
    in.paramstr = symb.paramstr;
	in.param = symb.param;
else
    in.paramstr = [];
    in.param = [];
end

%--------------------------------------------------------------------------
% objective function
%--------------------------------------------------------------------------
if isfield(symb,'Ob')

    linflagOb = false; % false only at the moment

    % calculate derivatives
	[obj,opts] = DTQP_IPFMINCON_symb(symb.Ob,in,linflagOb,false,opts);

    % initialize empty QP objective terms
    H = sparse([],[],[],in.nx,in.nx); f = sparse([],[],[],in.nx,1); c = 0;

    % check if there are OLQ elements
    if isfield(obj,'Ilin')

        % check if there are any quadratic objective function terms
        if ~isempty(obj.Ilin)

            % assign
            in.IDlin = obj.Ilin;
            in.IDnon = obj.Inon;

            % TODO: construct LQDO objective terms

        end
    end

    % assign
    in.obj = obj;

    % NLDO problem if any nonlinear elements
    % if ~isempty()
        lqdoflag = false;
    % end

else % only LQDO objective function terms

    % construct objective function matrices
    H = DTQP_createH(setup.L,setup.M,in,opts); % create Hessian
    f = DTQP_createf(setup.l,setup.m,in,opts); % create gradient
    c = DTQP_createc(setup.cL,setup.cM,in,opts); % determine constants

    % fmincon expects H directly
    H = H/2;

    % default field value
    in.obj = [];

end

%--------------------------------------------------------------------------
% dynamics
%--------------------------------------------------------------------------
if isfield(symb,'D')

    % determine nonlinear/linear dynamics
    [Aeq1,beq1,dyn,Idyn,in,opts] = DTQP_IPFMINCON_dyn(symb.D,in,linflag,opts,0,nt);

    % assign
    in.dyn = dyn;
    in.Ilambda.dyn = Idyn;

    % current number of nonlinear equality constraints
    nI = max(vertcat(in.Ilambda.dyn{:}));

    % NLDO problem if any nonlinear elements
    if ~isempty(Idyn)
        lqdoflag = false;
    end

else % only LQDO dynamic equations

    % construct linear defect constraint matrices
    [Aeq1,beq1,in] = DTQP_DEFECTS(setup.A,setup.B,setup.G,setup.d,in,opts);

    % default field value
    in.dyn = [];

end

%--------------------------------------------------------------------------
% general equality constraints
%--------------------------------------------------------------------------
if isfield(symb,'ceq')

    % determine nonlinear/linear equality constraints
    [Y,ceq,Iceq,opts] = DTQP_IPFMINCON_c(symb.ceq,in,linflag,opts,nI,nt);

    % combine
    setup.Y = [setup.Y;Y];

    % create linear path and boundary inequality constraints
    [Aeq2,beq2] = DTQP_create_YZ(setup.Y,in);

    % assign
    in.ceq = ceq;
    in.Ilambda.ceq = Iceq;

    % NLDO problem if any nonlinear elements
    if ~isempty(Iceq)
        lqdoflag = false;
    end

else % only LQDO equality constraints

    % create linear path and boundary equality constraints
    [Aeq2,beq2] = DTQP_create_YZ(setup.Y,in);

    % default field value
    in.ceq = [];

end

% combine linear equality constraints
Aeq = [Aeq1;Aeq2];
beq = [beq1;beq2]; % Aeq*X = beq

%--------------------------------------------------------------------------
% general inequality constraints
%--------------------------------------------------------------------------
if isfield(symb,'cin')

    % determine nonlinear/linear equality constraints
    [Z,cin,Icin,opts] = DTQP_IPFMINCON_c(symb.cin,in,linflag,opts,0,nt);

    % combine
    setup.Z = [setup.Z;Z];

    % create linear path and boundary inequality constraints
    [A,b] = DTQP_create_YZ(setup.Z,in); % A*X <= b

    % assign
    in.cin = cin;
    in.Ilambda.cin = Icin;

    % NLDO problem if any nonlinear elements
    if ~isempty(Icin)
        lqdoflag = false;
    end

else % only LQDO inequality constraints

    % create linear path and boundary inequality constraints
    [A,b] = DTQP_create_YZ(setup.Z,in);

    % default field value
    in.cin = [];

end

%--------------------------------------------------------------------------
% bounds
%--------------------------------------------------------------------------
[lb,ub] = DTQP_create_bnds(setup.LB,setup.UB,in);

% TODO: add flag
[A,b,Aeq,beq] = DTQP_scalingRows(A,b,Aeq,beq);

%--------------------------------------------------------------------------
% initial guess
%--------------------------------------------------------------------------
% TODO: create initial guess using DTQP_QLIN_guess.m
if isfield(in.p,'guess')

    % NOTE: the fields p.guess and p.Tguess need to be renamed

    % check different guess cases
    if size(in.p.guess,1) == 2

        % linear interpolation based on guess at end points
        X0 = interp1([in.t(1) in.t(end)],in.p.guess,in.t);

    else % arbitrary mesh provided

        % spline interpolation based on guess at end points
        X0 = interp1(in.p.Tguess,in.p.guess,in.t,'spline');

    end

    % extract
    X0uy = X0(:,1:(nu+ny));
    X0p = X0(1,end-np+1:end);

    % assign
    in.X0 = [X0uy(:);X0p(:)];

else % default guess
    in.X0 = ones(size(lb));
end

%--------------------------------------------------------------------------
% solve the optimization problem
%--------------------------------------------------------------------------
% check if this ended up as a LQDO problem
if lqdoflag

    % use quadprog instead of fmincon
    opts.solver.function = 'quadprog';

    % double Hessian for quadprog
    H = 2*H;

end

% solve
[X,F,in,opts] = DTQP_SOLVER(H,f,A,b,Aeq,beq,lb,ub,in,opts);

%--------------------------------------------------------------------------
% obtain outputs
%--------------------------------------------------------------------------
% return optimal controls, states, and parameters
T = in.t;
U = reshape(X((1:nu*nt)),nt,nu); % controls
Y = reshape(X((nu*nt+1:(nu+ny)*nt)),nt,ny); % states
P = reshape(X(((nu+ny)*nt+1:(nu+ny)*nt+np)),np,1); % parameters

% check for zero-order hold method and nan final controls
if strcmpi(opts.dt.defects,'ZO') || strcmpi(opts.dt.defects,'EF')
    U(end,:) = nan;
end

end

% dynamics
function [Aeq1,beq1,dyn,Idyn,in,opts] = DTQP_IPFMINCON_dyn(D,in,linflag,opts,nI,nt)

% calculate derivatives
[dyn,opts] = DTQP_IPFMINCON_symb(D,in,linflag,false,opts);

% number of constraints
nz = length(dyn.f);

% initialize
Idyn = cell(nz,1);

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

% additional constraints
function [YZ,c,Ic,opts] = DTQP_IPFMINCON_c(c,in,linflag,opts,nI,nt)

% check if pathboundary was provided
if isfield(c,'pathboundary')
    pathboundary = c.pathboundary;
end

% calculate derivatives
[c,opts] = DTQP_IPFMINCON_symb(c.func,in,linflag,true,opts);

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

% initialize
YZ = [];

% add to setup.Y structure
for k = 1:length(c.Ilin)

    % reset local index
    idx2 = 0;

    % states
    if ~isempty(c.A)

        % increment local index
        idx2 = idx2 + 1;

        % variable
        YZ(k).linear(idx2).right = 2; % states

        % assign
        YZ(k).linear(idx2).matrix = c.A(k,:);

    end

    % controls
    if ~isempty(c.B)

        % increment local index
        idx2 = idx2 + 1;

        % variable
        YZ(k).linear(idx2).right = 1; % controls

        % assign
        YZ(k).linear(idx2).matrix = c.B(k,:);

    end

    % parameters
    if ~isempty(c.G)

        % increment local index
        idx2 = idx2 + 1;

        % variable
        YZ(k).linear(idx2).right = 3; % parameters

        % assign
        YZ(k).linear(idx2).matrix = c.G(k,:);

    end

    % initial states
    if ~isempty(c.Ai)

        % increment local index
        idx2 = idx2 + 1;

        % variable
        YZ(k).linear(idx2).right = 4; % initial states

        % assign
        YZ(k).linear(idx2).matrix = c.Ai(k,:);

    end

    % final states
    if ~isempty(c.Af)

        % increment local index
        idx2 = idx2 + 1;

        % variable
        YZ(k).linear(idx2).right = 5; % final states

        % assign
        YZ(k).linear(idx2).matrix = c.Af(k,:);

    end

    % constant (negative disturbance)
    if ~isempty(c.d)

        % assign
        YZ(k).b = {'prod',-1,c.d(k,:)}; % constant

    end
end
end