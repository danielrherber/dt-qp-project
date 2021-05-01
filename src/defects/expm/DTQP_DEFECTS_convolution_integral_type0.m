%--------------------------------------------------------------------------
% DTQP_DEFECTS_convolution_integral_type0.m
% Compute the type 0 convolution integral:
% int(expm(A*(hk-s))*B(s+tk),s,0,hk)
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function I = DTQP_DEFECTS_convolution_integral_type0(A,B,in,~)

% extract
p = in.p; nt = in.nt; t = in.t; h = in.h;

% number of states
na = size(A,2);

% number of inputs
nb = size(B,2);

% initialize convolution output matrix
I = zeros(nt-1,na,nb);

% determine the type of matrix B
if isempty(B)
    return % no matrix provided
elseif iscell(B) || isa(B,'function_handle')
    B_type = 1; % time-varying matrix
else
    B_type = 0; % constant matrix
end

% check A type
if iscell(A) || isa(A,'function_handle') % time-varying matrix
    error('A should not be a function or time varying with the zero-order hold (ZOH) or first-order hold (FOH) methods')
elseif isempty(find(A,1))
    A_type = 0; % all zeros
elseif rank(full(A)) == na
    A_type = 1; % nonsingular constant A with some nonzeros
else
    A_type = 2; % singular constant A with some nonzeros
end

% tolerance for uniquetol (NOTE: potentially expose)
tol = 1e-12; % default

%--------------------------------------------------------------------------
% time-varying B matrix case
%--------------------------------------------------------------------------
if (B_type==1)

    % integration options
    options = {'ArrayValued',true,'RelTol',0,'AbsTol',1e-12};

    % go through each time segment/interval
    for k = 1:nt-1

        % extract segment/interval initial time and step size
        tk = t(k); hk = h(k);

        % check if any elements of A are nonzero
        if (A_type==0)
            integrand = @(s) shiftdim(DTQP_tmultiprod(B,p,s+tk));
        else % general case
            integrand = @(s) expm(A*(hk-s))*shiftdim(DTQP_tmultiprod(B,p,s+tk));
        end

        % compute integral
        I(k,:,:) = integral(integrand,0,hk,options{:});

    end

    % done
    return
end

%--------------------------------------------------------------------------
% nonsingular constant A with some nonzeros and constant B
% http://math.stackexchange.com/questions/658276/integral-of-matrix-exponential
%--------------------------------------------------------------------------
if (A_type==1)

    % find unique step size values
    [h_unique,~,IC_unique] = uniquetol(h,tol);

    % number of unique step sizes
    nh = length(h_unique);

    % initialize
    I_unique = zeros(nh,na,nb);

    % go through each unique step size
    for k = 1:nh

        % using method in the link above
        v = A\(expm(A*h_unique(k)) - eye(na))*B;

        % assign
        I_unique(k,:,:) = v;

    end

    % assign
    I = I_unique(IC_unique,:,:);

    % done
    return

end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% singular constant A with some nonzeros and constant B
%--------------------------------------------------------------------------
if (A_type==0) || (A_type==2)

    % find unique step size values
    [h_unique,~,IC_unique] = uniquetol(h,tol);

    % number of unique step sizes
    nh = length(h_unique);

    % initialize
    I_unique = zeros(nh,na,nb);

    % go through each unique step size
    for k = 1:nh

        % using c2d method
        v = expm([[A B]*h_unique(k); zeros(nb,na+nb)]);

        % assign
        I_unique(k,:,:) = v(1:na,na+1:na+nb);

    end

    % assign
    I = I_unique(IC_unique,:,:);

    % done
    return

end
%--------------------------------------------------------------------------

end