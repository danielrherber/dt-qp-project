%--------------------------------------------------------------------------
% DTQP_defects_TR_nonlin.m
% Create matrices for the trapezoidal method for nonlinear dynamics
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [Z,DZ] = DTQP_defects_TR_nonlin(X,dyn,in,opts,Dflag)

% extract
nu = in.nu; ny = in.ny; np = in.np; nt = in.nt; nX = in.nx;
p = in.p; t = in.t; h = in.h; ini = in.i; param = in.param;
f = dyn.f; Df = dyn.Df;

% states with nonlinear dynamics
if isfield(in,'IDnon')
	IDnon = in.IDnon;
else
    IDnon = 1:ny;
end

% nubmer of defect constraints
nz2 = length(IDnon);

% number of defect constraints per state
nz = nt - 1;

% initialize storage arrays
Isav = {}; Jsav = {}; Vsav = {};

% reshape optimization variables
P = X(end-np+1:end);
X = reshape(X(1:end-np),nt,[]);
P = repelem(P',nt,1);
X = [X,P];

% extract states for current optimization variables
Y = X(:,in.i{2});

% step size vector divided by 2
H2 = (h/2);

LR = repelem([1 2 3],[nu ny np]);
R = horzcat(ini{1:3});

%--------------------------------------------------------------------------
% compute defect constraints
%--------------------------------------------------------------------------
% calculate state derivative function values
fi = DTQP_qlin_update4tmatrix(f,[],X,param);
ft = DTQP_tmultiprod(fi,p,t);

% sum neighboring values
F = (ft + circshift(ft,[1 0]));

% extract relevant rows
F = F(2:end,:);

% calculate defect constraints
Z = -Y(2:end,IDnon) + Y(1:end-1,IDnon) + H2.*F;

% reshape to column vector
Z = Z(:);

%--------------------------------------------------------------------------
% compute Jacobian
%--------------------------------------------------------------------------
% check if Jacobian is requested
if ~Dflag
   DZ = [];
   return
end

% calculate Jacobian of state derivative function values
Dfi = DTQP_qlin_update4tmatrix(Df,[],X,param);
Dft = DTQP_tmultiprod(Dfi,p,t);

% go through each defect constraint
for ix = 1:nz2

    % go through each optimization variable (original form)
    for jx = R

        % extract
        V = Dft(:,ix,jx);
        
        % break if all entries are zero
        if ~any(V)
           continue 
        end
        
        % rows in DZ (defect constraints)
        I = ( (ix-1)*(nt-1) + 1:(ix*(nt-1)) )';

        % columns in DZ (optimization variables)
        J1 = DTQP_getQPIndex(R(jx),LR(jx),1,nt,nu,ny);

        % compute values
        V1 = H2.*V(1:end-1); % check H calculation
        V2 = H2.*V(2:end); % check H calculation

        % combine
        Isav{end+1} = I; Isav{end+1} = I;
        Jsav{end+1} = J1(1:end-1); Jsav{end+1} = J1(2:end);
        Vsav{end+1} = V1; Vsav{end+1} = V2;

    end

end

% combine
If = vertcat(Isav{:});
Jf = vertcat(Jsav{:});
Vf = vertcat(Vsav{:});

% create sparse matrix
FH = sparse(If,Jf,Vf,nz2*nz,nX);

% create initial matrix for state linear differences
ld = spdiags(repmat([1 -1],nz,1),[0 1],nz,nt);

% create nz copies block diagonal matrix
LDcell0 = repmat({sparse(0,nt)},ny,1);
LDcell0(IDnon) = {ld};
LD = blkdiag(LDcell0{:});

% combine with controls and parameters (empty matrices)
LD = [sparse(nz*nz2,nt*nu),LD,sparse(nz*nz2,np)];

% combine to create complete Jacobian for defect constraints
DZ = LD + FH;

end