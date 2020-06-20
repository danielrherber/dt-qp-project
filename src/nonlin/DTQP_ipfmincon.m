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
ldqoflag = false;
nI = 0;
in.Ilambda = [];
linflag = opts.qlin.olqflag;

% set default field value
if isfield(symb,'param')
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
	[obj,opts] = DTQP_IPFMINCON_symb(symb.Ob,in,linflagOb,opts);

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

    % calculate derivatives
	[dyn,opts] = DTQP_IPFMINCON_symb(symb.D,in,linflag,opts);

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
            end
            if ~isempty(dyn.B)
                setup.B = DTQP_QLIN_update_tmatrix(dyn.B,[],[],in.param);
            end
            if ~isempty(dyn.G)
                setup.G = DTQP_QLIN_update_tmatrix(dyn.G,[],[],in.param);
            end
            if ~isempty(dyn.d)
                setup.d = DTQP_QLIN_update_tmatrix(dyn.d,[],[],in.param);
            end

            % assign
            in.IDlin = dyn.Ilin;
            in.IDnon = dyn.Inon;

            % construct linear defect constraints
            [Aeq1,beq1,in] = DTQP_DEFECTS(setup.A,setup.B,setup.G,setup.d,in,opts);

        end
    end

    % assign
    in.dyn = dyn;
    in.Ilambda.dyn = Idyn;

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

    % calculate derivatives
	[ceq,opts] = DTQP_IPFMINCON_symb(symb.ceq,in,linflag,opts);

    pathboundary = true(length(ceq),1); % NOT CORRECT

    % number of constraints
    nz = length(ceq.f);

    % initialize
    Iceq = cell(nz,1);

    % go through each constraint
    for k = 1:nz

        % check if this is a path constraint
        if pathboundary
            Ii = nI + (1:nt)';
        else
            Ii = nI + 1;
        end

        % update
        nI = max(Ii);

        % assign
        Iceq{k} = Ii;

    end

    % create linear path and boundary inequality constraints
    [Aeq2,beq2] = DTQP_create_YZ(setup.Y,in);

    % assign
    in.ceq = ceq;
    in.Ilambda.ceq = Iceq;

else % only LQDO inequality constraints

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

    % calculate derivatives
	[cin,opts] = DTQP_IPFMINCON_symb(symb.cin,in,linflag,opts);

    pathboundary = true(length(cin),1); % NOT CORRECT

    % number of constraints
    nz = length(cin.f);

    % initialize
    Icin = cell(nz,1); nI = 0;

    % go through each constraint
    for k = 1:nz

        % check if this is a path constraint
        if pathboundary
            Ii = nI + (1:nt)';
        else
            Ii = nI + 1;
        end

        % update
        nI = max(Ii);

        % assign
        Icin{k} = Ii;

    end

    % create linear path and boundary inequality constraints
    [A,b] = DTQP_create_YZ(setup.Z,in); % A*X <= b

    % assign
    in.cin = cin;
    in.Ilambda.cin = Icin;

else % only LQDO inequality constraints

    % create linear path and boundary inequality constraints
    [A,b] = DTQP_create_YZ(setup.Z,in); % A*X <= b

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

    % interpolate initial guess
    X0 = interp1([in.t(1) in.t(end)],in.p.guess,in.t);

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
if ldqoflag

    % use quadprog instead of fmincon
    opts.qp.solver = 'quadprog';

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