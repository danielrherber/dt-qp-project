%--------------------------------------------------------------------------
% DTQP_IPFMINCON_hessian.m
% Compute Hessian of the Lagrangian
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function Ho = DTQP_IPFMINCON_hessian(X,lambda,obj,dyn,cin,ceq,Hin,in,opts)

% extract
nu = in.nu; ny = in.ny; np = in.np; ini = in.i; nx = in.nx;
p = in.p; t = in.t; nt = in.nt; param = in.param;
Ilambda = in.Ilambda; quadrature = opts.dt.quadrature;
derivativeflag = opts.method.derivatives; scaleflag = in.scaleflag;

% (potentially) apply linear scaling
if scaleflag

    % extract
    Xs = in.Xs; sm = in.sm; sc = in.sc;

    % unscale optimization variables
    Xunscaled = X.*sm + sc;

    % reshape optimization variables
    Xunscaled = DTQP_reshape_X(Xunscaled,np,nt,ini);

else

    % reshape optimization variables
    Xunscaled = DTQP_reshape_X(X,np,nt,ini);

end

% initialize row and column indices
LR = repelem([1 2 3 4 5],[nu ny np ny ny]);
R = horzcat(ini{1:5});
C = R;

% initialize storage arrays
Isav = {}; Jsav = {}; Vsav = {};

%--------------------------------------------------------------------------
% compute objective function Hessian terms
%--------------------------------------------------------------------------
if ~isempty(obj)

    % extract
    h = in.h; w = in.w;

    % integrate nonlinear term
    % TODO: add more methods
    switch upper(quadrature)
        case 'CEF'
            error(' ')
        case 'CTR'
            % augment step size vector
            h0 = [0;h/2];

            % compute circshifted step size sums
            H = (h0 + circshift(h0,[-1 0]));
        case 'CQHS'
            error(' ')
        case 'G'
            % compute product
            H = w*(in.tf - in.t0)/2;
        case 'CC'
            error(' ')
    end

    % calculate second derivative values
    D2ft = DTQP_hessian(obj,p,t,Xunscaled,param,derivativeflag,1);

    % go through each row entry in the original problem form
    for ix = 1:length(R)

        % go through each column entry in the original problem form
        for jx = 1:length(C)

            % get current values
            v = D2ft(:,ix,jx);

            % check there are nonzero entries
            if any(v)

                % compute step size vector and second derivative product
                v = v.*H;

                % Hessian row and column index sequences
                r = DTQP_getQPIndex(R(ix),LR(ix),1,nt,nu,ny);
                c = DTQP_getQPIndex(C(jx),LR(jx),1,nt,nu,ny);

                % main diagonal
                Isav{end+1} = r; % rows
                Jsav{end+1} = c; % columns
                Vsav{end+1} = v; % values

            end
        end
    end
end

%--------------------------------------------------------------------------
% compute defect constraints Hessian terms
%--------------------------------------------------------------------------
if ~isempty(dyn)

    % extract
    Ilambdadyn = Ilambda.dyn; h = in.h;

    % number of constraints
    nz = length(dyn.f);

    % extract and reshape multipliers
    lambda0 = full(lambda.eqnonlin(horzcat(Ilambdadyn{:})));

    % compute defect constraints
    % TODO: add more methods
    switch upper(opts.dt.defects)
        case 'ZO' % zero-order hold
            error(' ')
        case 'EF' % Euler forward
            error(' ')
        case 'TR' % trapezoidal
            % augment step size and multiplier vector                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  tors with initial zeros
            h0 = [0;h/2];
            lambda0 = vertcat(zeros(1,nz),lambda0);

            % compute multiplier and step size vector product
            h0lambda = lambda0.*h0;

            % compute circshifted multiplier sums
            h0lambdadyn = h0lambda + circshift(h0lambda,[-1,0]);
        case 'HS' % Hermite-Simpson
            error(' ')
        case 'RK4' % fourth-order Runge-Kutta
            error(' ')
        case 'PS' % pseudospectral (both LGL and CGL)
            h0lambdadyn = lambda0;
        case 'HUEN' % Heun's method
            error(' ')
        case 'MODEF' % Modified Euler method
            error(' ')
    end

    % go through each constraint
    for k = 1:nz

        % calculate second derivative values
        D2ft = DTQP_hessian(dyn,p,t,Xunscaled,param,derivativeflag,k);

        % continue if all derivatives are zero (so D2ft is empty)
        if isempty(D2ft)
           continue
        end

        % scale
        if scaleflag
            D2ft = D2ft./Xs(:,k);
        end

        % go through each row entry in the original problem form
        for ix = 1:length(R)

            % go through each column entry in the original problem form
            for jx = 1:length(C)

                % get current values
                v = D2ft(:,ix,jx);

                % check there are nonzero entries
                if any(v)

                    % compute multiplier and second derivative product
                    v = v.*h0lambdadyn(:,k);

                    % Hessian row and column index sequences
                    r = DTQP_getQPIndex(R(ix),LR(ix),1,nt,nu,ny);
                    c = DTQP_getQPIndex(C(jx),LR(jx),1,nt,nu,ny);

                    % main diagonal
                    Isav{end+1} = r; % rows
                    Jsav{end+1} = c; % columns
                    Vsav{end+1} = v; % values

                end
            end
        end
    end
end

%--------------------------------------------------------------------------
% compute general equality constraints Hessian terms
%--------------------------------------------------------------------------
if ~isempty(ceq)

    % extract
    Ilambda_ceq = Ilambda.ceq; pathboundary = ceq.pathboundary;

    % number of constraints
    nz = length(ceq.f);

    % extract multipliers
    lambda_eqnonlin = full(lambda.eqnonlin);

    % go through each constraint
    for k = 1:nz

        % calculate second derivative values
        D2ft = DTQP_hessian(ceq,p,t,Xunscaled,param,derivativeflag,k);

        % extract relevant multipliers
        lambda_ceq = lambda_eqnonlin(Ilambda_ceq{k});

        % go through each row entry in the original problem form
        for ix = 1:length(R)

            % go through each column entry in the original problem form
            for jx = 1:length(C)

                % get current values
                v = D2ft(:,ix,jx);

                % check there are nonzero entries
                if any(v)

                    % check if a path or boundary constraint
                    if pathboundary(k)

                        % compute multiplier and second derivative product
                        v = v.*lambda_ceq;

                        % Hessian row and column index sequences
                        r = DTQP_getQPIndex(R(ix),LR(ix),1,nt,nu,ny);
                        c = DTQP_getQPIndex(C(jx),LR(jx),1,nt,nu,ny);

                    else % boundary constraint

                        % compute multiplier and second derivative product
                        v = v(1).*lambda_ceq;

                        % Hessian row and column index sequences
                        r = DTQP_getQPIndex(R(ix),LR(ix),0,nt,nu,ny);
                        c = DTQP_getQPIndex(C(jx),LR(jx),0,nt,nu,ny);

                    end

                    % main diagonal
                    Isav{end+1} = r; % rows
                    Jsav{end+1} = c; % columns
                    Vsav{end+1} = v; % values

                end
            end
        end
    end
end

%--------------------------------------------------------------------------
% compute general inequality constraints Hessian terms
%--------------------------------------------------------------------------
if ~isempty(cin)

    % extract
    Ilambda_cin = Ilambda.cin; pathboundary = cin.pathboundary;

    % number of constraints
    nz = length(cin.f);

    % extract multipliers
    lambda_ineqnonlin = full(lambda.ineqnonlin);

    % go through each constraint
    for k = 1:nz

        % calculate second derivative values
        D2ft = DTQP_hessian(cin,p,t,Xunscaled,param,derivativeflag,k);

        % extract relevant multipliers
        lambda_cin = lambda_ineqnonlin(Ilambda_cin{k});

        % go through each row entry in the original problem form
        for ix = 1:length(R)

            % go through each column entry in the original problem form
            for jx = 1:length(C)

                % get current values
                v = D2ft(:,ix,jx);

                % check there are nonzero entries
                if any(v)

                    % check if a path or boundary constraint
                    if pathboundary(k)

                        % compute multiplier and second derivative product
                        v = v.*lambda_cin;

                        % Hessian row and column index sequences
                        r = DTQP_getQPIndex(R(ix),LR(ix),1,nt,nu,ny);
                        c = DTQP_getQPIndex(C(jx),LR(jx),1,nt,nu,ny);

                    else % boundary constraint

                        % compute multiplier and second derivative product
                        v = v(1).*lambda_cin;

                        % Hessian row and column index sequences
                        r = DTQP_getQPIndex(R(ix),LR(ix),0,nt,nu,ny);
                        c = DTQP_getQPIndex(C(jx),LR(jx),0,nt,nu,ny);

                    end

                    % main diagonal
                    Isav{end+1} = r; % rows
                    Jsav{end+1} = c; % columns
                    Vsav{end+1} = v; % values

                end
            end
        end
    end
end

%--------------------------------------------------------------------------
% combine
%--------------------------------------------------------------------------
% combine sequences
I = vertcat(Isav{:});
J = vertcat(Jsav{:});
V = vertcat(Vsav{:});

% construct sparse Hessian
Ho = sparse(I,J,V,nx,nx);

% scale
if scaleflag
    Ho = sm'.*Ho.*sm;
end

% add constant quadratic term
if ~isempty(find(Hin,1))
    Ho = Ho + 2*Hin;
end


end