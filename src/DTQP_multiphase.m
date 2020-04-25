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

% number of phases
nphs = length(setup);

% initialize storage arrays and constants
Hs = cell(nphs,1); fs = Hs; c = 0; SM = Hs; SC = Hs;
As = Hs; bs = Hs; Aeqs = Hs; beqs = Hs; lbs = Hs; ubs = Hs;
LALs = Hs; LARs = Hs; Lbs = Hs; LAeqLs = Hs; LAeqRs = Hs; Lbeqs = Hs;

% determine flags
scaleflag = isfield(setup,"scaling");
sqpflag = opts.qlin.sqpflag;
trustregionflag = opts.qlin.trustregionflag;
reorderflag = opts.qp.reorder;
multiphaseflag = nphs > 1;

%----------------------------------------------------------------------
% TASK: transcribe the problem for each phase and combine
%----------------------------------------------------------------------
% save current displevel
displevel = opts.general.displevel;
opts.general.displevel = 0;

% potentially start the timer
if (displevel > 0) % minimal
    tic % start timer
end

% go through each phase structure
for phs = 1:nphs
    % pass the options for the current phase
    phsopts = opts;
    phsopts.dt = opts.dt(phs);

    % transcribe the problem for the current phase
    [Hi,fi,ci,Ai,bi,Aeqi,beqi,lbi,ubi,setupi,in(phs),phsopts] = ...
        DTQP_create(setup(phs),phsopts);

    % transcribe linkage constraints (only in multiphase problems)
    if multiphaseflag
        % check if the linkage constraint fields are empty
        if ~isfield(setup,'LZ')
            setup(phs).LZ = [];
        end
        if ~isfield(setup,'LY')
            setup(phs).LY = [];
        end

        % transcribe the inequality linkage constraints for the current phase
        if phs == 1
            [LALi,Lbi,LARi,~] = DTQP_linkage(setup(phs).LZ,[],in(phs));
        else
            [LALi,Lbi,LARi,~] = DTQP_linkage(setup(phs).LZ,setup(phs-1).LZ,in(phs));
        end

        % transcribe the equality linkage constraints for the current phase
        if phs == 1
            [LAeqLi,Lbeqi,LAeqRi,~] = DTQP_linkage(setup(phs).LY,[],in(phs));
        else
            [LAeqLi,Lbeqi,LAeqRi,~] = DTQP_linkage(setup(phs).LY,setup(phs-1).LY,in(phs));
        end
    end

    % (optional) simple scaling
    if scaleflag
        if multiphaseflag
            error('Scaling is not currently supported for multiphase problems')
            [Hi,fi,ci,Ai,bi,Aeqi,beqi,lbi,ubi,LALi,LARi,Lbi,LAeqLi,LAeqRi,Lbeqi,in(phs),SM{phs},SC{phs}] = DTQP_scaling(...
                Hi,fi,ci,Ai,bi,Aeqi,beqi,lbi,ubi,LALi,LARi,Lbi,LAeqLi,LAeqRi,Lbeqi,in(phs),setupi.scaling);
        else
            [Hi,fi,ci,Ai,bi,Aeqi,beqi,lbi,ubi,~,~,~,~,~,~,in(phs),SM{phs},SC{phs}] = DTQP_scaling(...
                Hi,fi,ci,Ai,bi,Aeqi,beqi,lbi,ubi,[],[],[],[],[],[],in(phs),setupi.scaling);
        end
    end

    % (optional) add SQP penalty matrix
    if sqpflag
        if multiphaseflag
            error('SQP is not currently supported for multiphase problems')
        else
            if isfield(setup,'D2')
                % create sqp penalty matrix
                Hsqpi = DTQP_qlin_sqpMatrix(setup.D2,in,opts);

                % combine with original hessian
                if isempty(Hi)
                    Hi = Hsqpi;
                else
                    % combine
                    Hi = Hi + Hsqpi;

                    % flags (need to expose)
                    mirrorflag = true;
                    etol = sqrt(eps);

                    if mirrorflag
                        % check if the matrix is symmetric positive semidefinite
                        if eigs(Hi,1,'smallestreal') >= -etol
                            % disp('Matrix is symmetric positive definite')
                        else
                            % disp('Matrix is not symmetric positive definite')

                            % mirrored version
                            [Um,Tm] = schur(full(Hi));
                            Hi = Um*abs(Tm)*Um';
                        end

                        Hi = (Hi+Hi')/2; % make symmetric, then times 2 for 1/2*x'*H*x form
                    end
                end

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
        % store linear term
        if isempty(fi) % handle empty case
            fs{phs} = sparse([],[],[],in(phs).nx,1);
        else
            fs{phs} = fi;
        end

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

        Hs{phs} = Hi;
        As{phs} = Ai; bs{phs} = bi;
        Aeqs{phs} = Aeqi; beqs{phs} = beqi;
        LALs{phs} = LALi; LARs{phs} = LARi; Lbs{phs} = Lbi;
        LAeqLs{phs} = LAeqLi; LAeqRs{phs} = LAeqRi; Lbeqs{phs} = Lbeqi;
    else
        H = Hi; f = fi;
        A = Ai; b = bi;
        Aeq = Aeqi; beq = beqi;
        lb = lbi; ub = ubi;
    end

    % combine constants
    c = c + ci;

end

% combine problem elements (only for multiphase problems)
if multiphaseflag

    % combine to create complete equality linkage constraints
    LAeqL = blkdiag(LAeqLs{:});
    LAeqL = [LAeqL,sparse([],[],[],size(LAeqL,1),size(LAeqRs{end},2))];
    LAeqR = blkdiag(LAeqRs{:});
    LAeqR = [sparse([],[],[],size(LAeqR,1),size(LAeqLs{1},2)),LAeqR];
    LAeq = LAeqL + LAeqR;

    % combine to create complete inequality linkage constraints
    LAL = blkdiag(LALs{:});
    LAL = [LAL,sparse([],[],[],size(LAL,1),size(LARs{end},2))];
    LAR = blkdiag(LARs{:});
    LAR = [sparse([],[],[],size(LAR,1),size(LALs{1},2)),LAR];
    LA = LAL + LAR;

    % all problem elements
    H = blkdiag(Hs{:}); f = vertcat(fs{:});
    A = vertcat(blkdiag(As{:}),LA); b = vertcat(bs{:},Lbs{:});
    Aeq = vertcat(blkdiag(Aeqs{:}),LAeq); beq = vertcat(beqs{:},Lbeqs{:});
    lb = vertcat(lbs{:}); ub = vertcat(ubs{:});

end

% (optional) update problem elements based on trust region
if trustregionflag
    [beq,lb,ub,opts] = DTQP_SQP_trustregion(A,b,Aeq,beq,lb,ub,opts);
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
[X,F,in,opts] = DTQP_solver(H,f,A,b,Aeq,beq,lb,ub,in,opts);
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
    if strcmpi(opts.dt.defects,'ZO') || strcmpi(opts.dt.defects,'EF')
        U(end,:) = nan;
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