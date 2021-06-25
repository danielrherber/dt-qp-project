%--------------------------------------------------------------------------
% DTQP_multiphase.m
% Construct and solve a multiphase LQDO problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% TODO: add reordering (see DTQP_singlephase)
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [T,U,Y,P,F,in,opts] = DTQP_multiphase(setup,opts)

% handle multiple-interval pseudospectral method
if any(strcmpi({opts.dt.defects},'PS-MI'))
    [setup,opts] = DTQP_multi_interval_pseudospectral(setup,opts);
end

% number of phases
nphs = length(setup);

% initialize storage arrays and constants
Hs = cell(nphs,1); fs = Hs; c = 0; SM = Hs; SC = Hs;
As = Hs; bs = Hs; Aeqs = Hs; beqs = Hs; lbs = Hs; ubs = Hs;
LALs = Hs; LARs = Hs; LLbs = Hs; LRbs = Hs;
LAeqLs = Hs; LAeqRs = Hs; LLbeqs = Hs; LRbeqs = Hs;

% determine flags
scaleflag = isfield(setup,"scaling");
scalerowflag = opts.method.scalematrixrows;
reorderflag = opts.method.reordervariables;
multiphaseflag = nphs > 1;
if strcmpi(opts.method.form,'qlin')
    sqpflag = opts.method.sqpflag;
    trustregionflag = opts.method.trustregionflag;
else
    opts.method.sqpflag = false;
    sqpflag = false;
    opts.method.trustregionflag = false;
    trustregionflag = false;
end

%----------------------------------------------------------------------
% TASK: transcribe the problem for each phase and combine
%----------------------------------------------------------------------
% save current displevel
displevel = opts.general.displevel;
opts.general.displevel = 0;

% go through each phase structure
for phs = 1:nphs

    % pass the options for the current phase
    phsopts = opts;
    phsopts.dt = opts.dt(phs);

    % transcribe the problem for the current phase
    [Hi,fi,ci,Ai,bi,Aeqi,beqi,lbi,ubi,setup2(phs+1),in(phs),phsopts] = ...
        DTQP_create(setup(phs),phsopts);

    % NOTE: setup2 shifts phs by +1 compared to setup

    % transcribe linkage constraints (only in multiphase problems)
    if multiphaseflag

        % transcribe the inequality linkage constraints for the current phase
        [LLAi,LLbi,LRAi,LRbi] = DTQP_linkage(setup2(phs+1).LZ,setup2(phs).LZ,in(phs));

        % transcribe the equality linkage constraints for the current phase
        [LLAeqi,LLbeqi,LRAeqi,LRbeqi] = DTQP_linkage(setup2(phs+1).LY,setup2(phs).LY,in(phs));

    else

        % empty linkage elements
        LLAi = []; LLbi= []; LRAi= []; LRbi= [];
        LLAeqi= []; LLbeqi = []; LRAeqi = []; LRbeqi = [];

    end

    % (optional) linear transformation scaling
    if scaleflag
        [Hi,fi,ci,Ai,bi,Aeqi,beqi,lbi,ubi,LLAi,LRAi,LLbi,LRbi,LLAeqi,LRAeqi,LLbeqi,LRbeqi,in(phs),SM{phs},SC{phs}] = DTQP_scalingLinear(...
            Hi,fi,ci,Ai,bi,Aeqi,beqi,lbi,ubi,LLAi,LRAi,LLbi,LRbi,LLAeqi,LRAeqi,LLbeqi,LRbeqi,in(phs),setup2(phs+1).scaling);
    end

    % (optional) constraint row scaling
    if scalerowflag
        [Ai,bi,Aeqi,beqi] = DTQP_scalingRows(Ai,bi,Aeqi,beqi);
    end

    % (optional) add SQP penalty matrix
    if sqpflag
        if multiphaseflag
            error('SQP is not currently supported for multiphase problems')
        else
            if isfield(setup2(phs+1),'D2')

                % create sqp penalty matrix
                [Hi,Hsqpi] = DTQP_QLIN_sqp_matrix(Hi,setup2(phs+1).D2,in,opts);

            else
                Hsqpi = [];
            end
        end
    end

    % (optional) reordering of optimization variables
    if ~multiphaseflag
        if reorderflag
            [Hi,fi,ci,Ai,bi,Aeqi,beqi,lbi,ubi] = DTQP_reorder(in(phs),Hi,fi,ci,Ai,bi,Aeqi,beqi,lbi,ubi);
        end
    end

    % store problem elements
    if multiphaseflag

        % store lower bounds
        if isempty(lbi) % handle empty case
            lbs{phs} = -inf(in(phs).nx,1);
        else
            lbs{phs} = lbi;
        end

        % store upper bounds
        if isempty(ubi) % handle empty case
            ubs{phs} = inf(in(phs).nx,1);
        else
            ubs{phs} = ubi;
        end

        % assign
        Hs{phs} = Hi; fs{phs} = fi;
        As{phs} = Ai; bs{phs} = bi;
        Aeqs{phs} = Aeqi; beqs{phs} = beqi;
        LALs{phs} = LLAi; LARs{phs} = LRAi; LLbs{phs} = LLbi; LRbs{phs} = LRbi;
        LAeqLs{phs} = LLAeqi; LAeqRs{phs} = LRAeqi; LLbeqs{phs} = LLbeqi; LRbeqs{phs} = LRbeqi;

    else

        H = Hi; f = fi; A = Ai; b = bi;
        Aeq = Aeqi; beq = beqi; lb = lbi; ub = ubi;

    end

    % combine constants
    c = c + ci;

end

% combine problem elements (only for multiphase problems)
if multiphaseflag

    % combine to create complete equality linkage constraints
    LLAeq = blkdiag(LAeqLs{:});
    LLAeq = [LLAeq,sparse([],[],[],size(LLAeq,1),size(LAeqRs{end},2))];
    LRAeq = blkdiag(LAeqRs{:});
    LRAeq = [sparse([],[],[],size(LRAeq,1),size(LAeqLs{1},2)),LRAeq];
    LAeq = LLAeq + LRAeq;
    Lbeq = vertcat(LLbeqs{:}) + vertcat(LRbeqs{:});

    % combine to create complete inequality linkage constraints
    LAL = blkdiag(LALs{:});
    LAL = [LAL,sparse([],[],[],size(LAL,1),size(LARs{end},2))];
    LAR = blkdiag(LARs{:});
    LAR = [sparse([],[],[],size(LAR,1),size(LALs{1},2)),LAR];
    LA = LAL + LAR;
    Lb = vertcat(LLbs{:}) + vertcat(LRbs{:});

    % all problem elements
    H = blkdiag(Hs{:}); f = vertcat(fs{:});
    A = vertcat(blkdiag(As{:}),LA); b = vertcat(bs{:},Lb);
    Aeq = vertcat(blkdiag(Aeqs{:}),LAeq); beq = vertcat(beqs{:},Lbeq);
    lb = vertcat(lbs{:}); ub = vertcat(ubs{:});

end

% (optional) update problem elements based on trust region
if trustregionflag
    [beq,lb,ub,opts] = DTQP_SQP_trust_region(A,b,Aeq,beq,lb,ub,opts);
end

% previous displevel
opts.general.displevel = displevel;
%--------------------------------------------------------------------------
% END TASK: transcribe the problem for each phase and combine
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% TASK: solve the QP
%--------------------------------------------------------------------------
% solve the optimization problem
[X,F,in,opts] = DTQP_SOLVER(H,f,A,b,Aeq,beq,lb,ub,in,opts);
%--------------------------------------------------------------------------
% END TASK: solve the QP
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% TASK: obtain outputs
%--------------------------------------------------------------------------
if isempty(X)
    T = []; U = []; Y = []; P = [];
    return
end

% SQP penalty value
if sqpflag
    if ~isempty(Hsqpi)
        in.output.sqppenalty = X'*Hsqpi*(X/2);
    end
end

% (optional) restore ordering
if ~multiphaseflag % currently doesn't work for multiphase problems
    if reorderflag
        X = DTQP_reorder(in,X);
    end
end

% (optional) unscale solution
if scaleflag
    SM = vertcat(SM{:}); % scaling vector
    SC = vertcat(SC{:}); % scaling vector
    X = X.*SM + SC; % unscale optimization variables
end

% add the constant term to objective function
F = F + c;

% return optimal controls, states, and parameters
nx = 0;
T = cell(nphs,1); U = cell(nphs,1); Y = cell(nphs,1); P = cell(nphs,1);
for phs = 1:nphs

    % extract
    inphs = in(phs);
    nt = inphs.nt; nu = inphs.nu; ny = inphs.ny; np = inphs.np;

    % extract and reshape
    T{phs} = inphs.t;
    U{phs} = reshape(X(nx+(1:nu*nt)),nt,nu); % controls
    Y{phs} = reshape(X(nx+(nu*nt+1:(nu+ny)*nt)),nt,ny); % states
    P{phs} = reshape(X(nx+((nu+ny)*nt+1:(nu+ny)*nt+np)),np,1); % parameters

    % increment phase optimization variable number starting point
    nx = nx + inphs.nx;

end

% combine
T = vertcat(T{:}); U = vertcat(U{:}); Y = vertcat(Y{:}); P = vertcat(P{:});

% check for zero-order hold method and nan final controls
if nphs == 1
    if strcmpi(opts.dt(1).defects,'ZO') || strcmpi(opts.dt(1).defects,'EF')
        % U(end,:) = nan;
        U(end,:) = U(end-1,:); % use nearest control point
    end
end

% get unique time values (assumes continuity constraints are satisfied)
% [T, IT, ~] = unique(T);
%
% U = U(IT,:);
% Y = Y(IT,:);
%--------------------------------------------------------------------------
% END TASK: obtain outputs
%--------------------------------------------------------------------------

end