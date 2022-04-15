%--------------------------------------------------------------------------
% DTQP_NLP.m
% Prepare and solve NLDO problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [T,U,Y,P,F,in,opts] = DTQP_NLP(setup,opts)

% initialize some stuff
[setup,in] = DTQP_initialize(setup,opts.dt);

% check internal information is already available
% (likely from previous solve)
if isfield(setup,'internalinfo')
    in.internalinfo = setup.internalinfo; % copy
else
    in.internalinfo = []; % empty
end

% extract
nu = in.nu; ny = in.ny; np = in.np; nt = in.nt;
element = setup.element;

% initialize
lqdoflag = true;
in.Ilambda = [];
olqflag = opts.method.olqflag;

% set default field value
if isfield(element,'parameter_values')

    % convert time-varying cell to matrix on specified mesh
    % NOTE: this won't work for methods that need the parameters and
    % different points than the original mesh
    if isa(element.parameter_values,'cell')
        element.parameter_values = squeeze(DTQP_tmatrix(element.parameter_values,setup.auxdata,in.t));
    end
    in.parameter_list = element.parameter_list;
	in.param = element.parameter_values;
else
    in.parameter_list = [];
    in.param = [];
end

%--------------------------------------------------------------------------
% objective function
%--------------------------------------------------------------------------
% parse nonlinear/linear objective terms
[H,f,c,lqdoflag,in,opts] = DTQP_NLP_parse_objective(setup,in,opts);

%--------------------------------------------------------------------------
% dynamics
%--------------------------------------------------------------------------
% parse nonlinear/linear dynamics
[Aeq1,beq1,Idyn,in,opts,nI] = DTQP_NLP_parse_dynamics(setup,in,opts,olqflag,nt);

% NLDO problem if any nonlinear elements
if ~isempty(Idyn)
    lqdoflag = false;
end

%--------------------------------------------------------------------------
% general equality constraints
%--------------------------------------------------------------------------
% parse nonlinear/linear equality constraints
[Aeq2,beq2,Ih,in,opts,setup] = DTQP_NLP_parse_additional_constraints(setup,in,opts,olqflag,nt,nI,false);

% NLDO problem if any nonlinear elements
if ~isempty(Ih)
    lqdoflag = false;
end

%--------------------------------------------------------------------------
% combine linear equality constraints (Aeq*X = beq)
%--------------------------------------------------------------------------
Aeq = [Aeq1;Aeq2];
beq = [beq1;beq2];

%--------------------------------------------------------------------------
% general inequality constraints
%--------------------------------------------------------------------------
% parse nonlinear/linear equality constraints
[A,b,Ig,in,opts,setup] = DTQP_NLP_parse_additional_constraints(setup,in,opts,olqflag,nt,0,true);

% NLDO problem if any nonlinear elements
if ~isempty(Ig)
    lqdoflag = false;
end

%--------------------------------------------------------------------------
% bounds
%--------------------------------------------------------------------------
[lb,ub] = DTQP_create_bnds(setup.LB,setup.UB,in);

%--------------------------------------------------------------------------
% initial guess
%--------------------------------------------------------------------------
in = DTQP_guess(setup,in);

%--------------------------------------------------------------------------
% scaling
%--------------------------------------------------------------------------
% extract
s = setup.scaling;

% determine flags
scaleflag = ~isempty(s);
in.scaleflag = scaleflag;
scalerowflag = opts.method.scalematrixrows;

% (optional) apply linear transformation scaling
if scaleflag

    % apply scaling on olq elements
    [H,f,c,A,b,Aeq,beq,lb,ub,~,~,~,~,~,~,~,~,in,sm,sc] = ...
        DTQP_scalingLinear(H,f,c,A,b,Aeq,beq,lb,ub,[],[],[],[],[],[],[],[],in,s);

    % scale guess
    in.X0 = (in.X0 - sc)./sm;

    % obtain continuous variable scaling matrix
    SM = reshape(sm(1:end-np),nt,[]);

    % state variable locations
    Iy = in.i{2};

    % state variable scaling matrix for nonlinear state derivatives
    if isfield(in,'IDnon')
        Xs = SM(:,Iy(in.IDnon));
    else
        Xs = SM(:,Iy);
    end

    % assign
    in.Xs = Xs; in.sm = sm; in.sc = sc;

else

    % assign (for consistent structures)
    in.Xs = []; in.sm = []; in.sc = [];

end

% (optional) constraint row scaling
if scalerowflag
    [A,b,Aeq,beq] = DTQP_scalingRows(A,b,Aeq,beq);
end

%--------------------------------------------------------------------------
% check if this ended up as am LQDO problem
%--------------------------------------------------------------------------
if lqdoflag

    % use quadprog instead of fmincon
    opts.solver.function = 'quadprog';

    % double Hessian for quadprog
    H = 2*H;

end

%--------------------------------------------------------------------------
% solve the optimization problem
%--------------------------------------------------------------------------
[X,F,in,opts] = DTQP_SOLVER(H,f,A,b,Aeq,beq,lb,ub,in,opts);

%--------------------------------------------------------------------------
% obtain outputs
%--------------------------------------------------------------------------
% (optional) unscale optimization variables
if scaleflag
    X = X.*sm + sc;
end

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