%--------------------------------------------------------------------------
% DTQP_defects_EF.m
% Create matrices for the Euler forward method
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [Aeq,beq] = DTQP_defects_EF(A,B,G,d,in,opts)

% extract some of the variables
nu = in.nu; ny = in.ny; np = in.np; nd = in.nd; nx = in.nx;
p = in.p; nt = in.nt; t = in.t; h = in.h;

% matrix form of I in the formulas
K = kron(eye(ny),ones(nt-1,1));

% initialize storage arrays
Isav = {}; Jsav = {}; Vsav = {};

%--------------------------------------------------------------------------
% calculate matrices and sequencing vectors
%--------------------------------------------------------------------------
% find time dependent matrices
At = DTQP_tmultiprod(A,p,t);
Bt = DTQP_tmultiprod(B,p,t);
Gt = DTQP_tmultiprod(G,p,t);
dt = DTQP_tmultiprod(d,p,t);

Jy = DTQP_indexcolumns(nt,ny,nu); % optimization variable (column) locations
Jys = [Jy;Jy+1]; % combine to create paired optimization variable locations
Ty = Jy - nt*nu; % time indexing vectors
Hy = repmat(h,ny,1); % vector of time steps

if nu > 0
    Ju = DTQP_indexcolumns(nt,nu,0); % optimization variable (column) locations
    Tu = Ju; % time indexing vectors
    Hu = repmat(h,nu,1); % vector of time steps
end

if np > 0
    Jp = kron(nt*(nu+ny)+(1:np)', ones(nt-1,1)); % optimization variable (column) locations
    Tp = DTQP_indexcolumns(nt,np,0); % time indexing vectors
    Hp = repmat(h,np,1); % vector of time steps
end

if nd > 0
    Hd = h; % vector of time steps
    Td = 1:nt-1; % time indexing vectors
end
%--------------------------------------------------------------------------

% defect constraint of row continuous constraints
for i = 1:ny
    % current defect constraint row indices
    DefectIndices = reshape((i-1)*(nt-1)+1:i*(nt-1),[],1);

    %----------------------------------------------------------------------
    % controls
    %----------------------------------------------------------------------
    if nu > 0
        % defect constraint (row) locations
        Iu = repmat(DefectIndices,nu,1);

        % extract matrices
        Bv = reshape(Bt(:,i,:),[],1);

        % theta values
        V3 = -Hu.*Bv(Tu); % theta 3

        % combine
        Is = Iu; Js = Ju; Vs = V3;

        % remove zeros
        ZeroIndex = find(~Vs);
        Is(ZeroIndex) = []; Js(ZeroIndex) = []; Vs(ZeroIndex) = [];

        % combine
        Isav{end+1} = Is; Jsav{end+1} = Js; Vsav{end+1} = Vs;

    end
    %----------------------------------------------------------------------

    %----------------------------------------------------------------------
    % states
    %----------------------------------------------------------------------
    % if ny > 0 % there always is at least one state
        % defect constraint (row) locations
        Iy = repmat(DefectIndices,ny,1);

        % extract matrices
        Av = reshape(At(:,i,:),[],1);

        % theta values
        V1 = -K(:,i) - Hy.*Av(Ty); % theta 1
        V2 = K(:,i); % theta 2

        % combine
        Is = [Iy;Iy]; Js = Jys; Vs = [V1;V2];

        % remove zeros
        ZeroIndex = find(~Vs);
        Is(ZeroIndex) = []; Js(ZeroIndex) = []; Vs(ZeroIndex) = [];

        % combine
        Isav{end+1} = Is; Jsav{end+1} = Js; Vsav{end+1} = Vs;

    % end
    %----------------------------------------------------------------------

    %----------------------------------------------------------------------
    % parameters
    %----------------------------------------------------------------------
    if np > 0
        % defect constraint (row) locations
        Is = repmat(DefectIndices,np,1);

        % extract matrices
        Gv = reshape(Gt(:,i,:),[],1);

        % theta values
        Vs = -Hp.*Gv(Tp); % theta 5

        % combine
        Js = Jp;

        % remove zeros
        ZeroIndex = find(~Vs);
        Is(ZeroIndex) = []; Js(ZeroIndex) = []; Vs(ZeroIndex) = [];

        % combine
        Isav{end+1} = Is; Jsav{end+1} = Js; Vsav{end+1} = Vs;

    end
    %----------------------------------------------------------------------
end

% combine
If = vertcat(Isav{:});
Jf = vertcat(Jsav{:});
Vf = vertcat(Vsav{:});

% output sparse matrix
Aeq = sparse(If,Jf,Vf,ny*(nt-1),nx);

%--------------------------------------------------------------------------
% disturbance
%--------------------------------------------------------------------------
if nd > 0
    % initialize storage arrays
    Isav = {}; Vsav = {};

    % defect constraint of row continuous constraints
    for i = 1:ny
        % defect constraint (row) locations
        Is = reshape((i-1)*(nt-1)+1:i*(nt-1),[],1);

        % extract matrices
        dv = reshape(dt(:,i,:),[],1);

        % nu values
        Vs = Hd.*dv(Td); % nu

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
%--------------------------------------------------------------------------
end