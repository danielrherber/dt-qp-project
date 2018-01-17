%--------------------------------------------------------------------------
% DTQP_defects_ZO.m
% Create matrices for the zero-order hold method
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [Aeq,beq] = DTQP_defects_ZO(A,B,G,d,p,opts)

    % extract some of the variables in p
    nt = p.nt; nu = p.nu; ns = p.ns; np = p.np;
    nd = p.nd; h = p.h; nx = p.nx;
    
    % matrix form of I in the formulas
    K = kron(eye(ns),ones(nt-1,1));

    %------------------------------------------------------------------
    % calculate matrices
    %------------------------------------------------------------------
    % find time dependent matrices
    if isa(A,'double') % not time varying
        % initialize
        if strcmp(opts.NType,'ED')
            At = repmat(expm( A*h(1) ),1,1,nt-1);
        else
            At = zeros(ns,ns,nt-1);
            for i = 1:nt-1
                At(:,:,i) = expm( A*h(i) );
            end
        end
    else
       error('A matrix cannot be time varying with ZOH defect method') 
    end   
    Bt = DTQP_convolution(A,B,p,opts);
    Gt = DTQP_convolution(A,G,p,opts);
    dt = DTQP_convolution(A,d,p,opts);

    % permute
    At = permute(At,[1,3,2]);
    Bt = permute(Bt,[1,3,2]);
    Gt = permute(Gt,[1,3,2]);
    dt = permute(dt,[1,3,2]);
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

            % extract matrices
            Bv = reshape(Bt(i,:,:),[],1);

            % theta values
            V3 = -Bv; % theta 3

            % remove zeros
            ZeroIndex = (V3==0);
            I(ZeroIndex) = []; J(ZeroIndex) = []; V3(ZeroIndex) = [];

            % combine with 
            If = [If,I]; Jf = [Jf,J]; Vf = [Vf;V3];
 
        end
        %------------------------------------------------------------------
        
        %------------------------------------------------------------------
        % states
        %------------------------------------------------------------------
        % if ns > 0 % always is at least one state
            I = repmat(DefectIndices,1,ns); % current defect constraint row indices
            J = nu*nt+1:(nu+ns)*nt; % current optimization variable column indices
            J(nt:nt:end) = []; % remove endpoints
            
            % extract matrices
            Av = reshape(At(i,:,:),[],1);

            % theta values
            V1 = -Av; % theta 1
            V2 = K(:,i); % theta 2
            
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
            
            % extract matrices
            Gv = reshape(Gt(i,:,:),[],1);

            % theta values
            V = -Gv; % theta 5  

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

            % extract matrices
            dv = reshape(dt(i,:,:),[],1);

            % nu values
            V = dv; % nu

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