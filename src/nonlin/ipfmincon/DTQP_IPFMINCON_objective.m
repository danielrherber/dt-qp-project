%--------------------------------------------------------------------------
% DTQP_IPFMINCON_objective.m
% Compute objective function value and gradient
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [fo,go] = DTQP_IPFMINCON_objective(X,obj,in,opts,Hin,fin)

% extract
auxdata = in.auxdata; t = in.t; np = in.np; nt = in.nt; ini = in.i; param = in.param;
quadrature = opts.dt.quadrature; scaleflag = in.scaleflag;

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

%--------------------------------------------------------------------------
% compute objective
%--------------------------------------------------------------------------
% initialize objective function value
fo = 0;

% handle nonlinear term
if ~isempty(obj)

    % extract
    f = obj.f;

    % calculate objective function values
    fi = DTQP_QLIN_update_tmatrix(f,[],Xunscaled,param);
    ft = DTQP_tmultiprod(fi,auxdata,t);

    % integrate nonlinear term
    % TODO: add more methods
    switch upper(quadrature)
    case 'CEF'
        error(' ')
    case 'CTR'
        fo = fo + trapz(t,ft);
    case 'CQHS'
        error(' ')
    case 'G'
        fo = fo + (in.tf - in.t0)/2*in.w'*ft;
    case 'CC'
        error(' ')
    end

end

% handle quadratic term
if ~isempty(find(Hin,1))

    % compute intermediate value for the gradient
    XH = X'*Hin;

    % add to objective
    fo = fo + XH*X;

end

% handle linear term
if ~isempty(find(fin,1))

    % add to objective
	fo = fo + fin'*X;

end

%--------------------------------------------------------------------------
% compute gradient
%--------------------------------------------------------------------------
if nargout > 1

    % initial value (use linear term directly)
    go = fin';

    % handle nonlinear term
    if ~isempty(obj)

        % extract
        h = in.h; w = in.w;

        % calculate gradient of objective function values
        Dft = DTQP_jacobian(obj,auxdata,t,Xunscaled,param,opts.method.derivatives);

        % integrate nonlinear term
        % TODO: add more methods
        switch upper(quadrature)
            case 'CEF'
                error(' ')
            case 'CTR'
                % augment step size vector
                h0 = [0;h/2];

                % compute circshifted step size sums
                H0 = (h0 + circshift(h0,[-1 0]));

                % compute product
                Dft = Dft.*H0;

            case 'CQHS'
                error(' ')
            case 'G'
                % compute product
                Dft = Dft.*w*(in.tf - in.t0)/2;

            case 'CC'
                error(' ')
        end

        % extract submatrices
        Dft_UY = squeeze(Dft(:,:,1:(in.nu+in.ny))); % controls and states
        Dft_P = Dft(:,:,(in.nu+in.ny+1):(in.nu+in.ny+in.np)); % parameters

        % sum parameter jacobian
        Dft_P = sum(Dft_P,1);

        % combine
        Dft = [Dft_UY(:);Dft_P(:)]';

        % scale and add to gradient
        if scaleflag
            go = go + sm'.*Dft;
        else
            go = go + Dft;
        end

    end

    % handle quadratic term
    if ~isempty(find(Hin,1))

        % add to gradient
        go = go + 2*XH;

    end

end

end