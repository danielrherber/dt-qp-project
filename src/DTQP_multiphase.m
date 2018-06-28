%--------------------------------------------------------------------------
% DTQP_multiphase.m
% Construct and solve a multiphase LQDO problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% TODO: add scaling (see DTQP_singlephase)
% TODO: add reordering (see DTQP_singlephase)
%--------------------------------------------------------------------------
% Primary Contributor: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% Link: https://github.com/danielrherber/dt-qp-projec
%--------------------------------------------------------------------------
function [T,U,Y,P,F,in,opts] = DTQP_multiphase(setup,opts)

    % number of phases
    nphs = length(setup);
    
    % initialize storage arrays and constants
    Hs = cell(nphs,1); fs = Hs; c = 0;
    As = Hs; bs = Hs; Aeqs = Hs; beqs = Hs; lbs = Hs; ubs = Hs;
    LALs = Hs; LARs = Hs; Lbs = Hs;  LAeqLs = Hs; LAeqRs = Hs; Lbeqs = Hs;

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
        phsopts = rmfield(phsopts,'dt');
        phsopts.dt = opts.dt(phs);
        
        % transcribe the problem for the current phase
        [Hi,fi,ci,Ai,bi,Aeqi,beqi,lbi,ubi,~,in(phs),phsopts] = ...
            DTQP_create(setup(phs),phsopts);

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

        % store other terms
        Hs{phs} = Hi;
        As{phs} = Ai;
        bs{phs} = bi;
        Aeqs{phs} = Aeqi;
        beqs{phs} = beqi;

        % combine constants
        c = c + ci;
        
        %------------------------------------------------------------------
        % set left/right based on current phase
        if ~isfield(setup,'LZ')
            setup(phs).LZ = [];
        end
        left = setup(phs).LZ;
        if phs == 1
            right = []; % empty
        else
            right = setup(phs-1).LZ;
        end
        
        % transcribe the inequality linkage constraints for the current phase
        [LALi,Lbi,LARi,~] = DTQP_linkage(left,right,in(phs));
        
        % store
        LALs{phs} = LALi;
        LARs{phs} = LARi;
        Lbs{phs} = Lbi;
        %------------------------------------------------------------------
        
        %------------------------------------------------------------------
        % set left/right based on current phase
        if ~isfield(setup,'LY')
            setup(phs).LY = [];
        end
        left = setup(phs).LY;
        if phs == 1
            right = []; % empty
        else
            right = setup(phs-1).LY;
        end
        
        % transcribe the equality linkage constraints for the current phase
        [LAeqLi,Lbeqi,LAeqRi,~] = DTQP_linkage(left,right,in(phs));
        
        % store
        LAeqLs{phs} = LAeqLi;
        LAeqRs{phs} = LAeqRi;
        Lbeqs{phs} = Lbeqi;
        %------------------------------------------------------------------

    end

    % combine to create complete equality linkage constraints
    LAeqL = blkdiag(LAeqLs{:});
    LAeqL = [LAeqL,sparse([],[],[],size(LAeqL,1),size(LAeqRs{end},2))];
    LAeqR = blkdiag(LAeqRs{:});
    LAeqR = [sparse([],[],[],size(LAeqR,1),size(LAeqLs{1},2)),LAeqR];
    LAeq = LAeqL+LAeqR;
    
    % combine to create complete inequality linkage constraints
    LAL = blkdiag(LALs{:});
    LAL = [LAL,sparse([],[],[],size(LAL,1),size(LARs{end},2))];
    LAR = blkdiag(LARs{:});
    LAR = [sparse([],[],[],size(LAR,1),size(LALs{1},2)),LAR];
    LA = LAL+LAR;
    
    % combine
    H = blkdiag(Hs{:});
    f = vertcat(fs{:});
    A = vertcat(blkdiag(As{:}),LA);
    b = vertcat(bs{:},Lbs{:});
    Aeq = vertcat(blkdiag(Aeqs{:}),LAeq);
    beq = vertcat(beqs{:},Lbeqs{:});
    lb = vertcat(lbs{:});
    ub = vertcat(ubs{:});

    % end the timer
    if (displevel > 0) % minimal
        in(end).QPcreatetime = toc;
    end

    % display to the command window
    if (displevel > 1) % verbose
        disp(['QP creation time: ', num2str(in(end).QPcreatetime), ' s'])
    end

    % previous displevel
    opts.general.displevel = displevel;
    %----------------------------------------------------------------------
    % END TASK: transcribe the problem for each phase and combine
    %----------------------------------------------------------------------

    %----------------------------------------------------------------------
    % TASK: solve the QP
    %----------------------------------------------------------------------
    % solve the optimization problem
    [X,F,in,opts] = DTQP_solver(H,f,A,b,Aeq,beq,lb,ub,in,opts);
    %----------------------------------------------------------------------
    % END TASK: solve the QP
    %----------------------------------------------------------------------

    %----------------------------------------------------------------------
    % TASK: obtain outputs
    %----------------------------------------------------------------------
    % add the constant term to objective function
    F = F + c;

    % return optimal controls, states, and parameters
    nx = 0;
    T = []; U = []; Y = []; P = [];
    for phs = 1:length(setup)
        inphs = in(phs);
        T = [T; in(phs).t];
        U = [U; reshape(X(nx+(1:inphs.nu*inphs.nt)),inphs.nt,inphs.nu)]; % controls
        Y = [Y; reshape(X(nx+(inphs.nu*inphs.nt+1:(inphs.nu+inphs.ny)*inphs.nt)),inphs.nt,inphs.ny)]; % states
        P = [P; reshape(X(nx+((inphs.nu+inphs.ny)*inphs.nt+1:(inphs.nu+inphs.ny)*inphs.nt+inphs.np)),inphs.np,1)]; % parameters
        % increment phase optimization variable number starting point
        nx = nx + inphs.nx;
    end

    % get unique time values (assumes continuity constraints are satisfied)
    % [T, IT, ~] = unique(T);
    % 
    % U = U(IT,:);
    % Y = Y(IT,:);
    %----------------------------------------------------------------------
    % END TASK: obtain outputs
    %----------------------------------------------------------------------

end