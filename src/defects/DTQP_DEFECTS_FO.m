%--------------------------------------------------------------------------
% DTQP_DEFECTS_FO.m
% Create matrices for the first-order hold method
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [Aeq,beq] = DTQP_DEFECTS_FO(A,B,G,d,in,opts)

% extract some of the variables
nu = in.nu; ny = in.ny; np = in.np; nd = in.nd; nx = in.nx;
p = in.p; nt = in.nt; t = in.t; h = in.h;

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
At_expm = DTQP_DEFECTS_expm(A,in,opts);
Bt_type0 = DTQP_DEFECTS_convolution_integral_type0(A,B,in,opts);
Gt_type0 = DTQP_DEFECTS_convolution_integral_type0(A,G,in,opts);
dt_type0 = DTQP_DEFECTS_convolution_integral_type0(A,d,in,opts);
Bt_type1 = DTQP_DEFECTS_convolution_integral_type1(A,B,in,opts);

% check if time-dependent matrices are empty
Aflag = ~isempty(At_expm);
Bflag = ~isempty(Bt_type0) || ~isempty(Bt_type1);
Gflag = ~isempty(Gt_type0);
dflag = ~isempty(dt_type0);

% calculate matrix products
Jy = DTQP_DEFECTS_index_columns(nt,ny,nu); % optimization variable (column) locations
Jys = [Jy;Jy+1]; % combine to create paired optimization variable locations

if nu > 0
    Ju = DTQP_DEFECTS_index_columns(nt,nu,0); % optimization variable (column) locations
    Jus = [Ju;Ju+1]; % combine to create paired optimization variable locations
end

if np > 0
    Jp = kron(nt*(nu+ny)+(1:np)', ones(nt-1,1)); % optimization variable (column) locations
end

if nd > 0
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
        C2 = reshape(Bt_type0(:,i,:),[],1);
        C3 = reshape(Bt_type1(:,i,:),[],1);

        % check if any entries are nonzero
        if any(C2) || any(C3)

            % defect constraint (row) locations
            Iu = repmat(DefectIndices,nu,1);

            % theta values
            V3 = C3-C2; % theta 3
            V4 = -C3; % theta 4

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
            C1 = reshape(At_expm(:,i,:),[],1); % extract matrices
        else
            C1 = []; % empty
        end

        % check if any entries are nonzero
        if any(K(:,yi)) || any(C1)

            % defect constraint (row) locations
            Iy = repmat(DefectIndices,ny,1);

            % theta values
            V1 = -C1; % theta 1
            V2 = K(:,i); % theta 2

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
        Gv = reshape(Gt_type0(:,i,:),[],1);

        % check if any entries are nonzero
        if any(Gv)

            % defect constraint (row) locations
            Is = repmat(DefectIndices,np,1);

            % theta values
            Vs = -Gv; % theta 5

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
        dv = reshape(dt_type0(:,i,:),[],1);

        % check if any entries are nonzero
        if any(dv)

            % yi = IDlin(i);

            % defect constraint (row) locations
            Is = reshape((i-1)*(nt-1)+1:i*(nt-1),[],1);

            % nu values
            Vs = dv; % nu

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