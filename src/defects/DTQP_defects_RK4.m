%%--------------------------------------------------------------------------
% DTQP_defects_RK4.m
% Create matrices for the fourth-order Runge-Kutta method
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [Aeq,beq] = DTQP_defects_RK4(A,B,G,d,in,opts)

% extract some of the variables
nu = in.nu; ny = in.ny; np = in.np; nd = in.nd; nx = in.nx;
p = in.p; nt = in.nt; t = in.t; h = in.h; tm = in.tm;

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

% find matrix values for  time grid midpoints
Am = DTQP_tmultiprod(A,p,tm);
Bm = DTQP_tmultiprod(B,p,tm);
Gm = DTQP_tmultiprod(G,p,tm);
dm = DTQP_tmultiprod(d,p,tm);

% calculate matrix products
A0 = At(1:nt-1,:,:);
A1 = At(2:nt,:,:);
AmAm = multiprod(Am,Am,[2 3],[2 3]); % 1
AmAmA0 = multiprod(AmAm,A0,[2 3],[2 3]); % 2
A1AmAm = multiprod(A1,AmAm,[2 3],[2 3]); % 3
AmA0 = multiprod(Am,A0,[2 3],[2 3]); % 4
A1Am = multiprod(A1,Am,[2 3],[2 3]); % 5
A1AmAmA0 = multiprod(A1Am,AmA0,[2 3],[2 3]); % 6
Jy = DTQP_indexcolumns(nt,ny,nu); % optimization variable (column) locations
Jys = [Jy;Jy+1]; % combine to create paired optimization variable locations
Ty = Jy - nt*nu; % time indexing vectors
Hy = repmat(h,ny,1); % vector of time steps

if nu > 0
    B0 = Bt(1:nt-1,:,:);
    AmB0 = multiprod(Am,B0,[2 3],[2 3]); % 7
    AmBm = multiprod(Am,Bm,[2 3],[2 3]); % 8
    A1Bm = multiprod(A1,Bm,[2 3],[2 3]); % 9
    AmAmB0 = multiprod(AmAm,B0,[2 3],[2 3]); % 10
    A1AmBm = multiprod(A1Am,Bm,[2 3],[2 3]); % 11
    A1AmAmB0 = multiprod(A1,AmAmB0,[2 3],[2 3]); % 12
    Ju = DTQP_indexcolumns(nt,nu,0); % optimization variable (column) locations
    Jus = [Ju;Ju+1]; % combine to create paired optimization variable locations
    Tu = Ju; % time indexing vectors
    Hu = repmat(h,nu,1); % vector of time steps
end

if np > 0
    G0 = Gt(1:nt-1,:,:);
    AmG0 = multiprod(Am,G0,[2 3],[2 3]); % 13
    AmGm = multiprod(Am,Gm,[2 3],[2 3]); % 14
    A1Gm = multiprod(A1,Gm,[2 3],[2 3]); % 15
    AmAmG0 = multiprod(AmAm,G0,[2 3],[2 3]); % 16
    A1AmGm = multiprod(A1Am,Gm,[2 3],[2 3]); % 17
    A1AmAmG0 = multiprod(A1,AmAmG0,[2 3],[2 3]); % 18
    Jp = kron(nt*(nu+ny)+(1:np)', ones(nt-1,1)); % optimization variable (column) locations
    Tp = DTQP_indexcolumns(nt,np,0); % time indexing vectors
    Hp = repmat(h,np,1); % vector of time steps
end

if nd > 0
    d0 = dt(1:nt-1,:,:);
    Amd0 = multiprod(Am,d0,[2 3],[2 3]); % 19
    Amdm = multiprod(Am,dm,[2 3],[2 3]); % 20
    A1dm = multiprod(A1,dm,[2 3],[2 3]); % 21
    AmAmd0 = multiprod(AmAm,d0,[2 3],[2 3]); % 22
    A1Amdm = multiprod(A1Am,dm,[2 3],[2 3]); % 23
    A1AmAmd0 = multiprod(A1,AmAmd0,[2 3],[2 3]); % 24
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
        Bmv = reshape(Bm(:,i,:),[],1);
        AmB0v = reshape(AmB0(:,i,:),[],1);
        AmBmv = reshape(AmBm(:,i,:),[],1);
        A1Bmv = reshape(A1Bm(:,i,:),[],1);
        AmAmB0v = reshape(AmAmB0(:,i,:),[],1);
        A1AmBmv = reshape(A1AmBm(:,i,:),[],1);
        A1AmAmB0v = reshape(A1AmAmB0(:,i,:),[],1);

        % theta values
        V3 = -Hu.*( Bmv/3 + Bv(Tu)/6 + Hu.*AmB0v/6 + Hu.*AmBmv/12 + ...
            Hu.*A1Bmv/12 + Hu.^2.*AmAmB0v/12 + Hu.^2.*A1AmBmv/24 + ...
            Hu.^3.*A1AmAmB0v/24 ); % theta 3
        V4 = -Hu.*( Bmv/3 + Bv(Tu+1)/6 + Hu.*AmBmv/12 + Hu.*A1Bmv/12 + ...
            Hu.^2.*A1AmBmv/24 ); % theta 4

        % combine
        Is = [Iu;Iu]; Js = Jus; Vs = [V3;V4];

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
        Amv = reshape(Am(:,i,:),[],1);
        AmAmv = reshape(AmAm(:,i,:),[],1);
        AmAmA0v = reshape(AmAmA0(:,i,:),[],1);
        A1AmAmv = reshape(A1AmAm(:,i,:),[],1);
        AmA0v = reshape(AmA0(:,i,:),[],1);
        A1Amv = reshape(A1Am(:,i,:),[],1);
        A1AmAmA0v = reshape(A1AmAmA0(:,i,:),[],1);

        % theta values
        V1 = -K(:,i) - Hy.*( Amv*2/3 + Av(Ty+1)/6 + Av(Ty)/6 + ...
            Hy.*AmAmv/6 + Hy.*A1Amv/6 + Hy.*AmA0v/6 + Hy.^2.*A1AmAmv/12 + ...
            Hy.^2.*AmAmA0v/12 + Hy.^3.*A1AmAmA0v/24 ); % theta 1
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
        Gmv = reshape(Gm(:,i,:),[],1);
        AmG0v = reshape(AmG0(:,i,:),[],1);
        AmGmv = reshape(AmGm(:,i,:),[],1);
        A1Gmv = reshape(A1Gm(:,i,:),[],1);
        AmAmG0v = reshape(AmAmG0(:,i,:),[],1);
        A1AmGmv = reshape(A1AmGm(:,i,:),[],1);
        A1AmAmG0v = reshape(A1AmAmG0(:,i,:),[],1);

        % theta values
        V51 = -Hp.*( Gmv/3 + Gv(Tp)/6 + Hp.*AmG0v/6 + Hp.*AmGmv/12 + ...
            Hp.*A1Gmv/12 + Hp.^2.*AmAmG0v/12 + Hp.^2.*A1AmGmv/24 + ...
            Hp.^3.*A1AmAmG0v/24 );
        V52 = -Hp.*( Gmv/3 + Gv(Tp+1)/6 + Hp.*AmGmv/12 + Hp.*A1Gmv/12 + ...
            Hp.^2.*A1AmGmv/24 );
        Vs = V51 + V52; % theta 5

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
        dmv = reshape(dm(:,i,:),[],1);
        Amd0v = reshape(Amd0(:,i,:),[],1);
        Amdmv = reshape(Amdm(:,i,:),[],1);
        A1dmv = reshape(A1dm(:,i,:),[],1);
        AmAmd0v = reshape(AmAmd0(:,i,:),[],1);
        A1Amdmv = reshape(A1Amdm(:,i,:),[],1);
        A1AmAmd0v = reshape(A1AmAmd0(:,i,:),[],1);

        % nu values
        Vnu1 = Hd.*( dmv/3 + dv(Td)/6 + Hd.*Amd0v/6 + Hd.*Amdmv/12 + ...
            Hd.*A1dmv/12 + Hd.^2.*AmAmd0v/12 + Hd.^2.*A1Amdmv/24 + ...
            Hd.^3.*A1AmAmd0v/24 );
        Vnu2 = Hd.*( dmv/3 + dv(Td+1)/6 + Hd.*Amdmv/12 + Hd.*A1dmv/12 + ...
            Hd.^2.*A1Amdmv/24 );
        Vs = Vnu1 + Vnu2; % theta 5

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