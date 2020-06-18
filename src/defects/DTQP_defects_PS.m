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

% indices for linear defect constraints
if isfield(in,'IDlin')
	IDlin = in.IDlin;
else
    IDlin = 1:ny;
end

% number of linear defect constraints
nz = length(IDlin);

% differentiation matrix
D = in.D;
D = sparse(D);
Dns = cell(ny,1);
Dns(:) = {sparse(0,nt)};
Dns(IDlin) = {D};

% combine with controls and parameters (empty matrices)
Dadd = [sparse(nz*nt,nu*nt),blkdiag(Dns{:}),sparse(nz*nt,np)];

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

Jy = reshape(nu*nt+1:(nu+ny)*nt,[],1); % optimization variable (column) locations

if nu > 0
    Ju = reshape(1:nu*nt,[],1); % optimization variable (column) locations
end

if np > 0
    Jp = kron(nt*(nu+ny)+(1:np)',ones(nt,1)); % optimization variable (column) locations
end
%--------------------------------------------------------------------------

% defect constraint of row continuous constraints
for i = 1:nz

    % current defect constraint row indices
    DefectIndices = reshape((i-1)*nt+1:i*nt,[],1);

    %----------------------------------------------------------------------
    % controls
    %----------------------------------------------------------------------
    if nu > 0

        % extract matrices
        Vs = reshape(Bt(:,i,:),[],1);

        % check if any entries are nonzero
        if any(Vs)

            % defect constraint (row) locations
            Is = repmat(DefectIndices,nu,1);

            % combine
            Js = Ju;

            % combine
            Isav{end+1} = Is; Jsav{end+1} = Js; Vsav{end+1} = Vs;

        end
    end
    %----------------------------------------------------------------------

    %----------------------------------------------------------------------
    % states
    %----------------------------------------------------------------------
    if ny > 0

        % extract matrices
        Vs = reshape(At(:,i,:),[],1);

        % check if any entries are nonzero
        if any(Vs)

            % defect constraint (row) locations
            Is = repmat(DefectIndices,ny,1);

            % combine
            Js = Jy;

            % combine
            Isav{end+1} = Is; Jsav{end+1} = Js; Vsav{end+1} = Vs;

        end
    end
    %----------------------------------------------------------------------

    %----------------------------------------------------------------------
    % parameters
    %----------------------------------------------------------------------
    if np > 0

        % extract matrices
        Vs = reshape(Gt(:,i,:),[],1);

        % check if any entries are nonzero
        if any(Vs)

            % defect constraint (row) locations
            Is = repmat(DefectIndices,np,1);

            % combine
            Js = Jp;

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

% product with half time horizon length
Vf = -0.5*h*Vf;

% output sparse matrix
Aeq = sparse(If,Jf,Vf,nz*nt,nx) + Dadd;

%--------------------------------------------------------------------------
% disturbance
%--------------------------------------------------------------------------
if nd > 0

    % initialize storage arrays
    Isav = {}; Vsav = {};

    % defect constraint of row continuous constraints
    for i = 1:ny

        % extract matrices
        Vs = reshape(dt(:,i,:),[],1);

        % check if any entries are nonzero
        if any(Vs)

            % defect constraint (row) locations
            Is = reshape((i-1)*nt+1:i*nt,[],1);

            % combine
            Isav{end+1} = Is; Vsav{end+1} = Vs;

        end
    end

    % combine
    If = vertcat(Isav{:});
    Vf = vertcat(Vsav{:});

    % product with half time horizon length
    Vf = 0.5*h*Vf;

    % output sparse matrix
    beq = sparse(If,1,Vf,nz*nt,1);
else
    % output sparse matrix
    beq = sparse(nz*nt,1);
end
%--------------------------------------------------------------------------
end