%--------------------------------------------------------------------------
% DTQP_NLP_additional_constraints.m
% Compute nonlinear constraint values and Jacobian
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [G,DG] = DTQP_NLP_additional_constraints(X,con,in,opts,Dflag)

% extract
nu = in.nu; ny = in.ny; np = in.np; nt = in.nt; nX = in.nx;
auxdata = in.auxdata; t = in.t; ini = in.i; I_stored = in.I_stored; param = in.param;
f = con.f; pathboundary = con.pathboundary; scaleflag = in.scaleflag;

% (potentially) apply linear scaling
if scaleflag

    % extract
    sm = in.sm; sc = in.sc;

    % unscale optimization variables
    Xunscaled = X.*sm + sc;

    % reshape optimization variables
    Xunscaled = DTQP_reshape_X(Xunscaled,np,nt,ini);

else

    % reshape optimization variables
    Xunscaled = DTQP_reshape_X(X,np,nt,ini);

end

% number of constraints
nz = length(f);

%--------------------------------------------------------------------------
% compute constraint value
%--------------------------------------------------------------------------
% calculate constraint function values
fi = DTQP_QLIN_update_tmatrix(f,[],Xunscaled,param);
ft = DTQP_tmultiprod(fi,auxdata,t);

% initialize
G = cell(nz,1);

% go through each constraint
for ix = 1:nz

    % check if the constraint is a path or boundary constraint
    if pathboundary(ix)

        % path constraint
        G{ix} = ft(:,ix);

    else

        % boundary constraint
        G{ix} = ft(1,ix);

    end
end

% combine
G = vertcat(G{:});

%--------------------------------------------------------------------------
% compute Jacobian
%--------------------------------------------------------------------------
% check if Jacobian is requested
if ~Dflag
	DG = [];
	return
end

% initialize row and column indices
LR = repelem([1 2 3 4 5],[nu ny np ny ny]);
R = horzcat(ini{1:5});

% calculate Jacobian of the constraints
Dft = DTQP_jacobian(con,auxdata,t,Xunscaled,param,opts.method.derivatives);

% initialize storage arrays
Isav = {}; Jsav = {}; Vsav = {};

% initialize constraint row
r0 = 0;

% go through each constraint
for ix = 1:nz

    % initialize
    r = r0;

    % go through each column entry in the original problem form
    for jx = 1:length(R)

        % get current values
        v = Dft(:,ix,jx);

        % check there are nonzero entries
        if any(v)

            % check if the constraint is a path or boundary constraint
            if pathboundary(ix)

                % rows in DZ
                r = (r0+1:r0+nt)';

                % columns in DZ (optimization variables)
                c = DTQP_getQPIndex(R(jx),LR(jx),1,nt,I_stored);

            else % boundary constraint

                % row in DZ
                r = r0+1;

                % columns in DZ (optimization variables)
                c = DTQP_getQPIndex(R(jx),LR(jx),0,nt,I_stored);

                % only need first value
                v = v(1);

            end

            % main diagonal
            Isav{end+1} = r; % rows
            Jsav{end+1} = c; % columns
            Vsav{end+1} = v; % values

        end
    end

    % increment constraint row index
    r0 = r(end);

end

% combine
I = vertcat(Isav{:});
J = vertcat(Jsav{:});
V = vertcat(Vsav{:});

% create sparse matrix
DG = sparse(I,J,V,sum((pathboundary*nt)+(~pathboundary)),nX);

% scale
if scaleflag
    DG = sm'.*DG;
end

end