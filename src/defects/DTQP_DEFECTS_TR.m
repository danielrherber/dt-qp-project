%%--------------------------------------------------------------------------
% DTQP_DEFECTS_TR.m
% Create matrices for the trapezoidal method
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [Aeq,beq] = DTQP_DEFECTS_TR(A,B,G,d,in,opts)

% extract some of the variables
nu = in.nu; ny = in.ny; np = in.np; nd = in.nd; nx = in.nx;
auxdata = in.auxdata; nt = in.nt; t = in.t; h = in.h;

% states with linear dynamics
if isfield(in,'IDlin')
	IDlin = in.IDlin;
else
    IDlin = 1:ny;
end

% number of defect constraints
nz = length(IDlin);

% matrix form of I in the formulas
K = kron(eye(ny),ones(nt-1,1));

% initialize storage arrays
Isav = {}; Jsav = {}; Vsav = {};

%--------------------------------------------------------------------------
% calculate matrices and sequencing vectors
%--------------------------------------------------------------------------
% find time-dependent matrices
At = DTQP_tmultiprod(A,auxdata,t);
Bt = DTQP_tmultiprod(B,auxdata,t);
Gt = DTQP_tmultiprod(G,auxdata,t);
dt = DTQP_tmultiprod(d,auxdata,t);

% check if time-dependent matrices are empty
Aflag = ~isempty(At);
Bflag = ~isempty(Bt);
Gflag = ~isempty(Gt);
dflag = ~isempty(dt);

% calculate matrix products
Jy = DTQP_DEFECTS_index_columns(nt,ny,nu); % optimization variable (column) locations
Jys = [Jy;Jy+1]; % combine to create paired optimization variable locations
Ty = Jy - nt*nu; % time indexing vectors
Hy = repmat(0.5*h,ny,1); % vector of time steps

if nu > 0
    Ju = DTQP_DEFECTS_index_columns(nt,nu,0); % optimization variable (column) locations
    Jus = [Ju;Ju+1]; % combine to create paired optimization variable locations
    Tu = Ju; % time indexing vectors
    Hu = repmat(0.5*h,nu,1); % vector of time steps
end

if np > 0
    Jp = kron(nt*(nu+ny)+(1:np)', ones(nt-1,1)); % optimization variable (column) locations
    Tp = DTQP_DEFECTS_index_columns(nt,np,0); % time indexing vectors
    Hp = repmat(0.5*h,np,1); % vector of time steps
end

if nd > 0
    Hd = 0.5*h; % vector of time steps
    Td = 1:nt-1; % time indexing vectors
end
%--------------------------------------------------------------------------

% defect constraint of row continuous constraints
for i = 1:nz

    yi = IDlin(i);

    % current defect constraint row indices
    DefectIndices = reshape((i-1)*(nt-1)+1:i*(nt-1),[],1);

    %----------------------------------------------------------------------
    % controls
    %----------------------------------------------------------------------
    if nu > 0 && Bflag

        % extract matrices
        Bv = reshape(Bt(:,i,:),[],1);

        % check if any entries are nonzero
        if any(Bv)

            % defect constraint (row) locations
            Iu = repmat(DefectIndices,nu,1);

            % theta values
            V3 = -Hu.*Bv(Tu); % theta 3
            V4 = -Hu.*Bv(Tu+1); % theta 4

            % combine
            Is = [Iu;Iu]; Js = Jus; Vs = [V3;V4];

            % remove zeros
            ZeroIndex = ~Vs;
            Is(ZeroIndex) = []; Js(ZeroIndex) = []; Vs(ZeroIndex) = [];

            % combine
            Isav{end+1} = Is; Jsav{end+1} = Js; Vsav{end+1} = Vs;

        end
    end
    %----------------------------------------------------------------------

    %----------------------------------------------------------------------
    % states
    %----------------------------------------------------------------------
    if ny > 0

        % extract A matrix if nonempty
        if Aflag
            Av = reshape(At(:,i,:),[],1); % extract matrices
        else
            Av = []; % empty
        end

        % check if any entries are nonzero
        if any(K(:,yi)) || any(Av)

            % defect constraint (row) locations
            Iy = repmat(DefectIndices,ny,1);

            % theta values
            if any(Av)
                V1 = -K(:,yi) - Hy.*Av(Ty); % theta 1
                V2 = K(:,yi) - Hy.*Av(Ty+1); % theta 2
            else
                V1 = -K(:,yi); % theta 1
                V2 = K(:,yi); % theta 2
            end

            % combine
            Is = [Iy;Iy]; Js = Jys; Vs = [V1;V2];

            % remove zeros
            ZeroIndex = ~Vs;
            Is(ZeroIndex) = []; Js(ZeroIndex) = []; Vs(ZeroIndex) = [];

            % combine
            Isav{end+1} = Is; Jsav{end+1} = Js; Vsav{end+1} = Vs;

        end
    end
    %----------------------------------------------------------------------

    %----------------------------------------------------------------------
    % parameters
    %----------------------------------------------------------------------
    if np > 0 && Gflag

        % extract matrices
        Gv = reshape(Gt(:,i,:),[],1);

        % check if any entries are nonzero
        if any(Gv)

            % defect constraint (row) locations
            Is = repmat(DefectIndices,np,1);

            % theta values
            Vs = -Hp.*( Gv(Tp) + Gv(Tp+1) ); % theta 5

            % combine
            Js = Jp;

            % remove zeros
            ZeroIndex = ~Vs;
            Is(ZeroIndex) = []; Js(ZeroIndex) = []; Vs(ZeroIndex) = [];

            % combine
            Isav{end+1} = Is; Jsav{end+1} = Js; Vsav{end+1} = Vs;

        end
    end
    %----------------------------------------------------------------------
end

% combine
If = vertcat(Isav{:});
Jf = vertcat(Jsav{:});
Vf = vertcat(Vsav{:});

% output sparse matrix
Aeq = sparse(If,Jf,Vf,nz*(nt-1),nx);

%--------------------------------------------------------------------------
% disturbance
%--------------------------------------------------------------------------
if nd > 0  && dflag

    % initialize storage arrays
    Isav = {}; Vsav = {};

    % defect constraint of row continuous constraints
    for i = 1:nz

        % extract matrices
        dv = reshape(dt(:,i,:),[],1);

        % check if any entries are nonzero
        if any(dv)

            % yi = IDlin(i);

            % defect constraint (row) locations
            Is = reshape((i-1)*(nt-1)+1:i*(nt-1),[],1);

            % nu values
            Vs = Hd.*( dv(Td) + dv(Td+1) ); % nu

            % remove zeros
            ZeroIndex = ~Vs;
            Is(ZeroIndex) = []; Vs(ZeroIndex) = [];

            % combine
            Isav{end+1} = Is; Vsav{end+1} = Vs;
        end
    end

    % combine
    If = vertcat(Isav{:});
    Vf = vertcat(Vsav{:});

    % output sparse matrix
    beq = sparse(If,1,Vf,nz*(nt-1),1);
else
    % output sparse matrix
    beq = sparse(nz*(nt-1),1);
end
%--------------------------------------------------------------------------
end