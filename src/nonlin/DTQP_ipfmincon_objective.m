%--------------------------------------------------------------------------
% DTQP_ipfmincon_objective.m
% Compute objective function value and gradient
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [fo,go] = DTQP_ipfmincon_objective(X,obj,in,opts,Hin,fin)

% extract
p = in.p; t = in.t; np = in.np; nt = in.nt; param = in.param;
quadrature = opts.dt.quadrature;

% store initial optimization variable vector
Xo = X;

% reshape optimization variables
P = X(end-np+1:end);
X = reshape(X(1:end-np),nt,[]);
P = repelem(P',nt,1);
X = [X,P];

%--------------------------------------------------------------------------
% compute objective
%--------------------------------------------------------------------------
% initialize objective function value
fo = 0;

% handle nonlinear term
if isfield(obj,'f')

    % extract
    f = obj.f;

    % calculate state derivative function values
    fi = DTQP_qlin_update4tmatrix(f,[],X,param);
    ft = DTQP_tmultiprod(fi,p,t);

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
    XH = Xo'*Hin;

    % add to objective
    fo = fo + XH*Xo;

end

% handle linear term
if ~isempty(find(fin,1))

    % add to objective
	fo = fo + fin'*Xo;

end

%--------------------------------------------------------------------------
% compute gradient
%--------------------------------------------------------------------------
if nargout > 1

    % initial value (use linear term directly)
    go = fin';

    % handle nonlinear term
    if isfield(obj,'Df')

        % extract
        Df = obj.Df; h = in.h; w = in.w;

        % calculate state derivative function values
        Dfi = DTQP_qlin_update4tmatrix(Df,[],X,param);
        Dft = DTQP_tmultiprod(Dfi,p,t);

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

                % FIX:
                Dft(nt*(in.nu+in.ny)+np+1:end) = [];

            case 'CQHS'
                error(' ')
            case 'G'
                % compute product
                Dft = Dft.*w*(in.tf - in.t0)/2;

                % FIX:
                Dft(nt*(in.nu+in.ny)+np+1:end) = [];

            case 'CC'
                error(' ')
        end

        % ensure column vector
        Dft = Dft(:)';

        % add to gradient
        go = go + Dft;

    end

    % handle quadratic term
    if ~isempty(find(Hin,1))

        % add to gradient
        go = go + 2*XH;

    end

end
end