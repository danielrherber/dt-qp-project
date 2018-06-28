%--------------------------------------------------------------------------
% DTQP_defects_HS.m
% Create matrices for the Hermite-Simpson method
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [Aeq,beq] = DTQP_defects_HS(A,B,G,d,in,opts)

    % extract some of the variables
    nu = in.nu; ny = in.ny; np = in.np; nd = in.nd; nx = in.nx;
    p = in.p; nt = in.nt; t = in.t; h = in.h; tm = in.tm;

    % matrix form of I in the formulas
    K = kron(eye(ny),ones(nt-1,1));

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

    % find matrix values for time grid midpoints
    Am = DTQP_tmultiprod(A,p,tm);
    Bm = DTQP_tmultiprod(B,p,tm);
    Gm = DTQP_tmultiprod(G,p,tm);
    dm = DTQP_tmultiprod(d,p,tm);

    % time indexing vectors
    Im = 1:nt-1;
    Ip = 2:nt;

    % calculate matrix products
    AmA0 = multiprod(Am,At(Im,:,:),[2 3],[2 3]);
    AmA1 = multiprod(Am,At(Ip,:,:),[2 3],[2 3]);
    Jy = DTQP_indexcolumns(nt,ny,nu); % optimization variable (column) locations
    Jys = [Jy;Jy+1]; % combine to create paired optimization variable locations
    Ty = Jy - nt*nu; % time indexing vectors
    Hy = repmat(h,ny,1); % vector of time steps

    if nu > 0
        AmB0 = multiprod(Am,Bt(Im,:,:),[2 3],[2 3]);
        AmB1 = multiprod(Am,Bt(Ip,:,:),[2 3],[2 3]);
        Ju = DTQP_indexcolumns(nt,nu,0); % optimization variable (column) locations
        Jus = [Ju;Ju+1]; % combine to create paired optimization variable locations
        Tu = Ju; % time indexing vectors
        Hu = repmat(h,nu,1); % vector of time steps
    end

    if np > 0
        AbGm = multiprod(Am,Gt(Im,:,:),[2 3],[2 3]);
        AbGp = multiprod(Am,Gt(Ip,:,:),[2 3],[2 3]);
        Jp = kron(nt*(nu+ny)+(1:np)', ones(nt-1,1)); % optimization variable (column) locations
        Tp = DTQP_indexcolumns(nt,np,0); % time indexing vectors
        Hp = repmat(h,np,1); % vector of time steps
    end

    if nd > 0
        Abdm = multiprod(Am,dt(Im,:,:),[2 3],[2 3]);
        Abdp = multiprod(Am,dt(Ip,:,:),[2 3],[2 3]);
        Hd = h; % vector of time steps
        Td = 1:nt-1; % time indexing vectors
    end
    %----------------------------------------------------------------------

    % defect constraint of row continuous constraints
    for i = 1:ny
        % current defect constraint row indices
        DefectIndices = reshape((i-1)*(nt-1)+1:i*(nt-1),[],1);

        %------------------------------------------------------------------
        % controls
        %------------------------------------------------------------------
        if nu > 0
            % defect constraint (row) locations
            Iu = repmat(DefectIndices,nu,1);

            % extract matrices
            Bv = reshape(Bt(:,i,:),[],1);
            Bmv = reshape(Bm(:,i,:),[],1);
            AmBmv = reshape(AmB0(:,i,:),[],1);
            AmBpv = reshape(AmB1(:,i,:),[],1);

            % theta values           
            V3 = -Hu.*(Bv(Tu)/6 + Bmv/3 + Hu.*AmBmv/12); % theta 3
            V4 = -Hu.*(Bv(Tu+1)/6 + Bmv/3 - Hu.*AmBpv/12); % theta 4

            % combine
            Is = [Iu;Iu]; Js = Jus; Vs = [V3;V4];

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
            Iy = repmat(DefectIndices,ny,1);

            % extract matrices
            Av = reshape(At(:,i,:),[],1);
            Amv = reshape(Am(:,i,:),[],1);
            AmAmv = reshape(AmA0(:,i,:),[],1);
            AmApv = reshape(AmA1(:,i,:),[],1);

            % theta values
            V1 = -K(:,i) - Hy.*( Av(Ty)/6 + Amv/3 + Hy.*AmAmv/12 ); % theta 1
            V2 = K(:,i) - Hy.*( Av(Ty+1)/6 + Amv/3 - Hy.*AmApv/12 ); % theta 2

            % combine
            Is = [Iy;Iy]; Js = Jys; Vs = [V1;V2];

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
            Gmv = reshape(Gm(:,i,:),[],1);
            AmGmv = reshape(AbGm(:,i,:),[],1);
            AmGpv = reshape(AbGp(:,i,:),[],1);

            % theta values
            Vs = -Hp.*( Gv(Tp)/6 + Gmv*2/3 + Gv(Tp+1)/6 + Hp.*AmGmv/12 - Hp.*AmGpv/12 ); % theta 5

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
    Aeq = sparse(If,Jf,Vf,ny*(nt-1),nx);

    %------------------------------------------------------------------
	% disturbance
    %------------------------------------------------------------------
	if nd > 0
        % initialize storage arrays
        Isav = {}; Vsav = {};

        % defect constraint of row continuous constraints
        for i = 1:ny
            % defect constraint (row) locations
            Is = reshape((i-1)*(nt-1)+1:i*(nt-1),[],1);

            % extract matrices
            dv = reshape(dt(:,i,:),[],1);
            dmv = reshape(dm(:,i,:),[],1);
            Amdmv = reshape(Abdm(:,i,:),[],1);
            Amdpv = reshape(Abdp(:,i,:),[],1);

            % nu values
            Vs = Hd.*( dv(Td)/6 + 2*dmv/3 + dv(Td+1)/6 + Hd.*Amdmv/12 - Hd.*Amdpv/12 ); % nu

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
        beq = sparse(If,1,Vf,ny*(nt-1),1);   
    else
        % output sparse matrix
        beq = sparse([],[],[],ny*(nt-1),1);
	end
    %------------------------------------------------------------------
end