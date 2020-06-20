%--------------------------------------------------------------------------
% DTQP_IPFMINCON_additional_constraints.m
% Compute nonlinear constraint values and Jacobian
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [G,DG] = DTQP_IPFMINCON_additional_constraints(X,con,in,opts,Dflag)

% extract
nu = in.nu; ny = in.ny; np = in.np; nt = in.nt; nX = in.nx;
p = in.p; t = in.t; ini = in.i; param = in.param;
f = con.f; Df = con.Df;

% reshape optimization variables
P = X(end-np+1:end);
X = reshape(X(1:end-np),nt,[]);
P = repelem(P',nt,1);
X = [X,P];

% initialize row and column indices
LR = repelem([1 2 3],[nu ny np]);
R = horzcat(ini{1:3});

%--------------------------------------------------------------------------
% compute constraint value
%--------------------------------------------------------------------------
% calculate constraint function values
fi = DTQP_QLIN_update_tmatrix(f,[],X,param);
ft = DTQP_tmultiprod(fi,p,t);

% ensure column vector
G = ft(:);

%--------------------------------------------------------------------------
% compute Jacobian
%--------------------------------------------------------------------------
% check if Jacobian is requested
if ~Dflag
   DG = [];
   return
end

% calculate Jacobian of the constraints
Dfi = DTQP_QLIN_update_tmatrix(Df,[],X,param);
Dft = DTQP_tmultiprod(Dfi,p,t);

% number of constraints
nz = length(f);

% initialize storage arrays
Isav = {}; Jsav = {}; Vsav = {};

% go through each constraint
for ix = 1:nz

    % go through each column entry in the original problem form
    for jx = R

        % get current values
        v = Dft(:,ix,jx);

        % check there are nonzero entries
        if any(v)

            % rows in DZ (defect constraints)
            r = ( (ix-1)*nt + 1:(ix)*nt )';

            % columns in DZ (optimization variables)
            c = DTQP_getQPIndex(R(jx),LR(jx),1,nt,nu,ny);

            % main diagonal
            Isav{end+1} = r; % rows
            Jsav{end+1} = c; % columns
            Vsav{end+1} = v; % values

        end
    end
end

% combine
I = vertcat(Isav{:});
J = vertcat(Jsav{:});
V = vertcat(Vsav{:});

% create sparse matrix
DG = sparse(I,J,V,nz*nt,nX);

end