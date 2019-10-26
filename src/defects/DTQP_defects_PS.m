%--------------------------------------------------------------------------
% DTQP_defects_PS.m
% Create matrices for the pseudospectral method
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [Aeq,beq] = DTQP_defects_PS(A,B,G,d,in,opts)

    % extract some of the variables
    nu = in.nu; ny = in.ny; np = in.np; nd = in.nd; nx = in.nx;
    p = in.p; nt = in.nt; t = in.t; h = in.tf - in.t0;

    % differentiation matrix
    D = in.D;
    D = sparse(D);
    Dns = cell(1,ny);
    Dns(:) = {D};
    Dadd = [sparse([],[],[],ny*nt,nu*nt),blkdiag(Dns{:}),sparse([],[],[],ny*nt,np)];

    % initialize storage arrays
    Isav = {}; Jsav = {}; Vsav = {};

    %----------------------------------------------------------------------
    % calculate matrices and sequencing vectors
    %----------------------------------------------------------------------
    % find time dependent matrices
    At = DTQP_tmultiprod(A,p,t);
    Bt = DTQP_tmultiprod(B,p,t);
    Gt = DTQP_tmultiprod(G,p,t);
    dt = DTQP_tmultiprod(d,p,t);
    
    Jy = reshape(nu*nt+1:(nu+ny)*nt,[],1); % optimization variable (column) locations

    if nu > 0
        Ju = reshape(1:nu*nt,[],1); % optimization variable (column) locations
    end

    if np > 0
        Jp = kron(nt*(nu+ny)+(1:np)', ones(nt,1)); % optimization variable (column) locations
    end
    %----------------------------------------------------------------------

    % defect constraint of row continuous constraints
    for i = 1:ny
        % current defect constraint row indices
        DefectIndices = reshape((i-1)*nt+1:i*nt,[],1);

        %------------------------------------------------------------------
        % controls
        %------------------------------------------------------------------
        if nu > 0
            % defect constraint (row) locations
            Is = repmat(DefectIndices,nu,1);

            % extract matrices
            Bv = reshape(Bt(:,i,:),[],1);

            % values
            Vs = -0.5*h*Bv;

            % combine
            Js = Ju;

            % remove zeros
            ZeroIndex = find(~Vs);
            Is(ZeroIndex) = []; Js(ZeroIndex) = []; Vs(ZeroIndex) = [];

            % combine 
            Isav{end+1} = Is; Jsav{end+1} = Js; Vsav{end+1} = Vs;

        end
        %------------------------------------------------------------------

        %------------------------------------------------------------------
        % states
        %------------------------------------------------------------------
        % if ny > 0 % there always is at least one state
            % defect constraint (row) locations
            Is = repmat(DefectIndices,ny,1);

            % extract matrices
            Av = reshape(At(:,i,:),[],1);

            % values
            Vs = -0.5*h*Av;

            % combine
            Js = Jy;

            % remove zeros
            ZeroIndex = find(~Vs);
            Is(ZeroIndex) = []; Js(ZeroIndex) = []; Vs(ZeroIndex) = [];

            % combine 
            Isav{end+1} = Is; Jsav{end+1} = Js; Vsav{end+1} = Vs;

        % end
        %------------------------------------------------------------------

        %------------------------------------------------------------------
        % parameters
        %------------------------------------------------------------------
        if np > 0
            % defect constraint (row) locations
            Is = repmat(DefectIndices,np,1);

            % extract matrices
            Gv = reshape(Gt(:,i,:),[],1);

            % values
            Vs = -0.5*h*Gv;

            % combine
            Js = Jp;

            % remove zeros
            ZeroIndex = find(~Vs);
            Is(ZeroIndex) = []; Js(ZeroIndex) = []; Vs(ZeroIndex) = [];

            % combine 
            Isav{end+1} = Is; Jsav{end+1} = Js; Vsav{end+1} = Vs;

        end
        %------------------------------------------------------------------
    end

    % combine
    If = vertcat(Isav{:});
    Jf = vertcat(Jsav{:});
    Vf = vertcat(Vsav{:});

	% output sparse matrix   
    Aeq = sparse(If,Jf,Vf,ny*nt,nx) + Dadd;

    %------------------------------------------------------------------
	% disturbance
    %------------------------------------------------------------------
	if nd > 0
        % initialize storage arrays
        Isav = {}; Vsav = {};

        % defect constraint of row continuous constraints
        for i = 1:ny
            % defect constraint (row) locations
            Is = reshape((i-1)*nt+1:i*nt,[],1);

            % extract matrices
            dv = reshape(dt(:,i,:),[],1);

            % values
            Vs = 0.5*h*dv;

            % remove zeros
            ZeroIndex = find(~Vs);
            Is(ZeroIndex) = []; Vs(ZeroIndex) = [];

            % combine 
            Isav{end+1} = Is; Vsav{end+1} = Vs;

        end

        % combine
        If = vertcat(Isav{:});
        Vf = vertcat(Vsav{:});

        % output sparse matrix
        beq = sparse(If,1,Vf,ny*nt,1);   
    else
        % output sparse matrix
        beq = sparse([],[],[],ny*nt,1);
	end
    %------------------------------------------------------------------
end