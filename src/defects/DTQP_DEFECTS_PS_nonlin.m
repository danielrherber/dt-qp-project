%--------------------------------------------------------------------------
% DTQP_DEFECTS_PS_nonlin.m
% Create matrices for the pseudospectral method for nonlinear dynamics
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [Z,DZ] = DTQP_DEFECTS_PS_nonlin(X,dyn,in,opts,Dflag)

% extract
nu = in.nu; ny = in.ny; np = in.np; nt = in.nt; nX = in.nx;
auxdata = in.auxdata; t = in.t; D = in.D; ini = in.i; I_stored = in.I_stored;
f = dyn.f; scaleflag = in.scaleflag; param = in.param;

% (potentially) apply linear scaling
if scaleflag

    % extract
    Xs = in.Xs; sm = in.sm; sc = in.sc;

    % unscale optimization variables
    Xunscaled = X.*sm + sc;

    % reshape optimization variables
    X = DTQP_reshape_X(X,np,nt,ini);
    Xunscaled = DTQP_reshape_X(Xunscaled,np,nt,ini);

else

    % reshape optimization variables
    X = DTQP_reshape_X(X,np,nt,ini);

    % reshape optimization variables
    Xunscaled = X;

end

% compute time horizon length
h = in.tf - in.t0;

% states with nonlinear dynamics
if isfield(in,'IDnon')
	IDnon = in.IDnon;
else
    IDnon = 1:ny;
end

% number of defect constraints
nz2 = length(IDnon);

% number of defect constraints per state
nz = nt;

% initialize storage arrays
Isav = {}; Jsav = {}; Vsav = {};

% extract states for current optimization variables
Y = X(:,in.i{2});

LR = repelem([1 2 3],[nu ny np]);
R = horzcat(ini{1:3});

%--------------------------------------------------------------------------
% compute defect constraints
%--------------------------------------------------------------------------
% calculate state derivative function values
fi = DTQP_QLIN_update_tmatrix(f,[],Xunscaled,param);
ft = DTQP_tmultiprod(fi,auxdata,t);

% scale
if scaleflag
    ft = ft./Xs;
end

% calculate defect constraints
Z = ft - 2/h*D*Y(:,IDnon);

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
Dft = DTQP_jacobian(dyn,auxdata,t,Xunscaled,param,opts.method.derivatives);

% scale
if scaleflag
    Dft = Dft./Xs;
end

% go through each defect constraint
for ix = 1:nz2

    % go through each optimization variable (original form)
    for jx = R

        % compute values
        V = Dft(:,ix,jx);

        % break if all values are zero
        if ~any(V)
            continue
        end

        % rows in DZ (defect constraints)
        I = ( (ix-1)*nz + 1:(ix*nz) )' ;

        % columns in DZ (optimization variables)
        J = DTQP_getQPIndex(R(jx),LR(jx),1,nt,I_stored); % Hessian row index sequence

        % combine
        Isav{end+1} = I; Jsav{end+1} = J; Vsav{end+1} = V;

    end
end

% combine
If = vertcat(Isav{:});
Jf = vertcat(Jsav{:});
Vf = vertcat(Vsav{:});

% create sparse matrix
FH = sparse(If,Jf,Vf,nz2*nz,nX);

% scale
if scaleflag
    FH = sm'.*FH;
end

% differentiation matrix
D = sparse(-2/h*D);
Dns = repmat({sparse(0,nt)},ny,1);
Dns(IDnon) = {D};

% combine with controls and parameters (empty matrices)
LD = [sparse(nz2*nt,nu*nt),blkdiag(Dns{:}),sparse(nz2*nt,np)];

% combine to create complete Jacobian for defect constraints
DZ = FH + LD;

end