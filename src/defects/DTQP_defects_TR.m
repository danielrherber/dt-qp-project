%--------------------------------------------------------------------------
% DTQP_defects_TR.m
% Create matrices for the trapezoidal method
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [Aeq,beq] = DTQP_defects_TR(A,B,G,d,p,opts)

    % extract some of the variables in p
    nt = p.nt; nu = p.nu; ns = p.ns; np = p.np;
    nd = p.nd; h = p.h; nx = p.nx;

    % matrix form of I in the formulas
    K = kron(eye(ns),ones(nt-1,1));

    %------------------------------------------------------------------
    % calculate matrices
    %------------------------------------------------------------------
    % find time dependent matrices
    At = DTQP_tmultiprod(A,p);
    Bt = DTQP_tmultiprod(B,p);
    Gt = DTQP_tmultiprod(G,p);
    dt = DTQP_tmultiprod(d,p);
    %------------------------------------------------------------------

    % initialize sequences 
    If = []; Jf = []; Vf = [];

    % defect constraint of row continuous constraints
    for i = 1:ns
        % current defect constraint row indices
        DefectIndices = (i-1)*(nt-1)+1:i*(nt-1);

        %------------------------------------------------------------------
        % controls
        %------------------------------------------------------------------
        if nu > 0
            I = repmat(DefectIndices,1,nu); % current defect constraint row indices
            J = 1:nu*nt; % current optimization variable column indices
            J(nt:nt:nu*nt) = []; % remove endpoints
            % T = 1:nu*nt; T(nt:nt:nu*nt) = [];
            T = J; % time indexing vector
            H = repmat(h,nu,1); % vector of time steps

            % extract matrices
            Bv = reshape(Bt(:,i,:),[],1);

            % theta values
            V3 = -H.*Bv(T)/2; % theta 3
            V4 = -H.*Bv(T+1)/2; % theta 4

            % combine with
            Is = [I,I];
            Js = [J,J+1];
            Vs = [V3;V4];

            % remove zeros
            ZeroIndex = (Vs==0);
            Is(ZeroIndex) = []; Js(ZeroIndex) = []; Vs(ZeroIndex) = [];

            % combine with 
            If = [If,Is]; Jf = [Jf,Js]; Vf = [Vf;Vs];
 
        end
        %------------------------------------------------------------------

        %------------------------------------------------------------------
        % states
        %------------------------------------------------------------------
        % if ns > 0 % always is at least one state
            I = repmat(DefectIndices,1,ns); % current defect constraint row indices
            J = nu*nt+1:(nu+ns)*nt; % current optimization variable column indices
            J(nt:nt:end) = []; % remove endpoints
            % T = 1:ns*nt; T(nt:nt:ns*nt) = [];
            T = J - nt*nu; % time indexing vector (faster than line above)
            H = repmat(h,ns,1); % vector of time steps

            % extract matrices
            Av = reshape(At(:,i,:),[],1);

            % theta values
            V1 = -K(:,i) - H.*Av(T)/2; % theta 1
            V2 = K(:,i) - H.*Av(T+1)/2; % theta 2

            % combine
            Is = [I,I];
            Js = [J,J+1];
            Vs = [V1;V2];

            % remove zeros
            ZeroIndex = (Vs==0);
            Is(ZeroIndex) = []; Js(ZeroIndex) = []; Vs(ZeroIndex) = [];

            % combine 
            If = [If,Is]; Jf = [Jf,Js]; Vf = [Vf;Vs];
        % end
        %------------------------------------------------------------------

        %------------------------------------------------------------------
        % parameters
        %------------------------------------------------------------------
        if np > 0
            I = repmat(DefectIndices,1,np); % current defect constraint row indices
            J = kron(nt*(nu+ns)+(1:np), ones(1,nt-1)); % current optimization variable column indices
            T = 1:np*nt; T(nt:nt:np*nt) = []; % time indexing vector
            H = repmat(h,np,1); % vector of time steps

            % extract matrices
            Gv = reshape(Gt(:,i,:),[],1);

            % theta values
            V = -H/2.*( Gv(T) + Gv(T+1) ); % theta 5  

            % remove zeros
            ZeroIndex = (V==0);
            I(ZeroIndex) = []; J(ZeroIndex) = []; V(ZeroIndex) = [];

            % combine
            If = [If,I]; Jf = [Jf,J]; Vf = [Vf;V];
        end
        %------------------------------------------------------------------
    end

	% output sparse matrix   
    Aeq = sparse(If,Jf,Vf,ns*(nt-1),nx);

    %------------------------------------------------------------------
	% disturbance
    %------------------------------------------------------------------
    if nd > 0
        % initialize sequences 
        Ifb = []; Vfb = [];

        for i = 1:ns % defect constraint of row continuous constraints
            I = (i-1)*(nt-1)+1:i*(nt-1); % row (continuous)
            T = 1:nt-1; % time indexing vector
            H = h; % vector of time steps

            % extract matrices
            dv = reshape(dt(:,i,:),[],1);

            % nu values
            V = H/2.*( dv(T) + dv(T+1) ); % nu

            % remove zeros
            ZeroIndex = (V==0);
            I(ZeroIndex) = []; V(ZeroIndex) = [];

            % combine with 
            Ifb = [Ifb,I]; Vfb = [Vfb;V];

        end

        beq = sparse(Ifb,1,Vfb,ns*(nt-1),1);   
    else
        % output sparse matrix
        beq = sparse([],[],[],ns*(nt-1),1);
    end
    %------------------------------------------------------------------
end