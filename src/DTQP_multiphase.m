%--------------------------------------------------------------------------
% DTQP_multiphase.m
% Construct and solve a multiphase LQDO problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% TODO: parameters should be consistent across the phases
% TODO: add general linear phase linking (continuity) constraints
% TODO: add scaling (see DTQP_solve)
% TODO: add reordering (see DTQP_solve)
%--------------------------------------------------------------------------
% Primary Contributor: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% Link: https://github.com/danielrherber/dt-qp-projec
%--------------------------------------------------------------------------
function [T,U,Y,P,F,p,opts] = DTQP_multiphase(SETUP,opts)

    % initialize some stuff
    [SETUP,opts] = DTQP_default_opts(SETUP,opts);

    % initialize QP matrices
    H = []; f = []; c = 0; A = []; b = []; 
    Aeq = []; beq = []; lb = []; ub = [];
    
    %----------------------------------------------------------------------
    % TASK: transcribe the problem for each phase and combine
    %----------------------------------------------------------------------
    % save current displevel
    displevel = opts.displevel;
    opts.displevel = 0;

    % potentially start the timer
    if (displevel > 0) % minimal
        tic % start timer
    end
    
    % go through each phase structure
    for idx = 1:length(SETUP)
        % transcribe the problem
        [Hi,fi,ci,Ai,bi,Aeqi,beqi,lbi,ubi,~,p(idx),opts] =...
            DTQP_create(SETUP(idx),opts);

        % combine matrices
        H = [H, sparse([],[],[],size(H,1),size(Hi,2));...
            sparse([],[],[],size(Hi,1),size(H,2)), Hi];
        Aeq = [Aeq, sparse([],[],[],size(Aeq,1),size(Aeqi,2));...
            sparse([],[],[],size(Aeqi,1),size(Aeq,2)), Aeqi];
        A = [A, sparse([],[],[],size(A,1),size(Ai,2));...
            sparse([],[],[],size(Ai,1),size(A,2)), Ai];
        
        % combine vectors
        if isempty(fi) % handle empty case
            f = [f; sparse([],[],[],p(idx).nx,1)];
        else
            f = [f; fi];
        end
        
        if isempty(lbi) && isempty(ubi) % handle empty case
            lb = [lb; -inf(p(idx).nx,1)];
            ub = [ub; inf(p(idx).nx,1)];
        else
            lb = [lb; lbi];
            ub = [ub; ubi];
        end
        
        b = [b; bi];
        beq = [beq; beqi];

        % combine constants
        c = c + ci;
    end

    % add continuity constraints
    nx = 0;
    for idx = 1:length(SETUP)-1
        % left phase state indices (end)
        pL = p(idx);
        IL = nx + (pL.nu*pL.nt+pL.nt:pL.nt:(pL.nu+pL.ns)*pL.nt);
        nx = nx + pL.nx;
        
        % right phase state indices (initial)
        pR = p(idx+1);
        IR = nx + (pR.nu*pR.nt+1:pR.nt:(pR.nu+pR.ns)*pR.nt);
        
        % 
        Aeqi = sparse([1:length(IL),1:length(IR)],[IL,IR],...
            [ones(1,length(IL)),-ones(1,length(IR))],length(IR),size(Aeq,2));
        
        % combine
        Aeq = [Aeq; Aeqi];
        beq = [beq; sparse([],[],[],size(Aeqi,1),1)];
    end

    % end the timer
    if (displevel > 0) % minimal
        opts.QPcreatetime = toc;
    end

    % display to the command window
    if (displevel > 1) % verbose
        disp(['QP creation time: ', num2str(opts.QPcreatetime), ' s'])
    end

    % previous displevel
    opts.displevel = displevel;
    %----------------------------------------------------------------------
    % END TASK: transcribe the problem for each phase and combine
    %----------------------------------------------------------------------

    %----------------------------------------------------------------------
    % TASK: solve the QP
    %----------------------------------------------------------------------
    % solve the optimization problem
    [X,F,opts] = DTQP_solver(H,f,A,b,Aeq,beq,lb,ub,opts);
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
    for idx = 1:length(SETUP)
        pI = p(idx);
        T = [T;p(idx).t];
        U = [U; reshape(X(nx+(1:pI.nu*pI.nt)),pI.nt,pI.nu)]; % controls
        Y = [Y; reshape(X(nx+(pI.nu*pI.nt+1:(pI.nu+pI.ns)*pI.nt)),pI.nt,pI.ns)]; % states
        P = [P; reshape(X(nx+((pI.nu+pI.ns)*pI.nt+1:(pI.nu+pI.ns)*pI.nt+pI.np)),pI.np,1)]; % parameters
        % increment phase optimization variable number starting point
        nx = nx + pI.nx;
    end

    % get unique time values (assumes continuity constraints are satisfied)
    [T, IT, ~] = unique(T);

    U = U(IT,:);
    Y = Y(IT,:);
    %----------------------------------------------------------------------
    % END TASK: obtain outputs
    %----------------------------------------------------------------------

end