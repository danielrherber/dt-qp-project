%--------------------------------------------------------------------------
% DTQP_defects_PS.m
% Create matrices for the pseudospectral method
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [Aeq,beq] = DTQP_defects_PS(A,B,G,d,p,opts)

    % extract some of the variables in p
    nt = p.nt; nu = p.nu; ns = p.ns; np = p.np;
    nd = p.nd; h = p.tf - p.t0; nx = p.nx;

    % differentiation matrix
    D = p.D;
    D = sparse(D);
    Dns = cell(1,ns);
    Dns(:) = {D};
    Dadd = [sparse([],[],[],ns*nt,nu*nt),blkdiag(Dns{:}),sparse([],[],[],ns*nt,np)];

    %------------------------------------------------------------------
    % calculate matrices
    %------------------------------------------------------------------
    % find time dependent matrices
    At = DTQP_tmultiprod(A,p);
    Bt = DTQP_tmultiprod(B,p);
    Gt = DTQP_tmultiprod(G,p);
    dt = DTQP_tmultiprod(d,p);
    %------------------------------------------------------------------

    % determine locations and matrix values at this points
    If = []; Jf = [];  Vf = [];

    for i = 1:ns % defect constraint of row continuous constraints
        % current defect constraint row indices
        DefectIdx = (i-1)*nt+1:i*nt;

        %------------------------------------------------------------------
        % controls
        %------------------------------------------------------------------
        if nu ~= 0
            I = repmat(DefectIdx,1,nu); % row (continuous)
            J = 1:nu*nt; % column (continuous)
            Bv = reshape(Bt(:,i,:),[],1);
            V = -h/2*Bv;

            % remove zeros
            ZeroIndex = (V==0);
            I(ZeroIndex) = [];
            J(ZeroIndex) = [];
            V(ZeroIndex) = [];

            % combine with 
            If = [If,I];
            Jf = [Jf,J];
            Vf = [Vf;V];
        end
        %------------------------------------------------------------------

        %------------------------------------------------------------------
        % states
        %------------------------------------------------------------------
        % if ns ~= 0 % always is at least one state
            I = repmat(DefectIdx,1,ns); % row (continuous)
            J = nu*nt+1:(nu+ns)*nt; % column (continuous)
            Av = reshape(At(:,i,:),[],1);
            V = -h/2*Av;

            % remove zeros
            ZeroIndex = (V==0);
            I(ZeroIndex) = [];
            J(ZeroIndex) = [];
            V(ZeroIndex) = [];

            % combine with
            If = [If,I];
            Jf = [Jf,J];
            Vf = [Vf;V];
        % end
        %------------------------------------------------------------------

        %------------------------------------------------------------------
        % parameters
        %------------------------------------------------------------------
        if np ~= 0
            I = repmat(DefectIdx,1,ns); % row (continuous)
            J = kron(nt*(nu+ns)+(1:np), ones(1,nt));
            Gv = reshape(Gt(:,i,:),[],1);
            V = -h/2*Gv;

            % remove zeros
            ZeroIndex = (V==0);
            I(ZeroIndex) = [];
            J(ZeroIndex) = [];
            V(ZeroIndex) = [];

            % combine with 
            If = [If,I];
            Jf = [Jf,J];
            Vf = [Vf;V];
        end
        %------------------------------------------------------------------
    end

    % output sparse matrix
    Aeq = sparse(If,Jf,Vf,ns*nt,nx) + Dadd;

    %------------------------------------------------------------------
	% disturbance
    %------------------------------------------------------------------
    if nd ~= 0
        % initialize sequences 
        Ifb = []; Vfb = [];

        for i = 1:ns % defect constraint of row continuous constraints
            I = (i-1)*nt+1:i*nt; % row (continuous)
            dv = reshape(dt(:,i,:),[],1);
            V = h/2*dv;

            % remove zeros
            ZeroIndex = (V==0);
            I(ZeroIndex) = [];
            V(ZeroIndex) = [];

            % combine with 
            Ifb = [Ifb,I];
            Vfb = [Vfb;V];

        end

        beq = sparse(Ifb,1,Vfb,ns*nt,1);   
    else
        % output sparse matrix
        beq = sparse([],[],[],ns*nt,1);
    end
    %------------------------------------------------------------------

end