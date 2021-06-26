%--------------------------------------------------------------------------
% DTQP_DEFECTS_convolution_integral_type1.m
% Compute the type 1 convolution integral:
% int(s*expm(A*(hk-s))*B(s+tk),s,0,hk)/hk
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function I = DTQP_DEFECTS_convolution_integral_type1(A,B,in,~)

% extract
auxdata = in.auxdata; nt = in.nt; t = in.t; h = in.h;

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

% integration options
options = {'ArrayValued',true,'RelTol',0,'AbsTol',1e-12};

% tolerance for uniquetol (NOTE: potentially expose)
tol = 1e-12; % default

%--------------------------------------------------------------------------
% time-varying B matrix case
%--------------------------------------------------------------------------
if (B_type==1)

    % go through each time segment/interval
    for k = 1:nt-1

        % extract segment/interval initial time and step size
        tk = t(k); hk = h(k);

        % check if any elements of A are nonzero
        if (A_type==0)
            integrand = @(s) s*eye(na)*shiftdim(DTQP_tmultiprod(B,auxdata,s+tk));
        else % general case
            integrand = @(s) s*expm(A*(hk-s))*shiftdim(DTQP_tmultiprod(B,auxdata,s+tk));
        end

        % compute integral
        I(k,:,:) = integral(integrand,0,hk,options{:})/hk;

    end

    % done
    return

end

%--------------------------------------------------------------------------
% nonsingular constant A with some nonzeros and constant B
%--------------------------------------------------------------------------
if (A_type==1)

    % find unique step size values
    [h_unique,IA_unique,IC_unique] = uniquetol(h,tol);

    % number of unique step sizes
    nh = length(h_unique);

    % initialize
    I_unique = zeros(nh,na,nb);

    % compute simplier integral (with nonsingular constant A and B matrices)
    I_minus_t = shiftdim(DTQP_DEFECTS_convolution_integral_type0(A,B,in,[]),1);

    % go through each unique step size
    for k = 1:nh

        % extract segment/interval step size
        hk = h(k);

        % extract
        if ndims(I_minus_t) == 3
            int_QPp = I_minus_t(:,:,IA_unique(k));
        else
            int_QPp = I_minus_t(:,IA_unique(k));
        end

        % formula based on matrix integration by parts (see bottom)
        I_unique(k,:,:) = A\(int_QPp/hk - B);

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

        % extract segment/interval initial time and step size
        hk = h_unique(k);

        % check if any elements of A are nonzero
        if (A_type==0)
            % I_unique(k,:,:) = integral(@(s) s*eye(na)*B,0,hk,options{:})/hk;
            I_unique(k,:,:) = hk/2*eye(na)*B; % closed-form
        else % general case
            I_unique(k,:,:) = integral(@(s) s*expm(A*(hk-s))*B,0,hk,options{:})/hk;
        end

    end

    % assign
    I = I_unique(IC_unique,:,:);

    % done
    return

end
%--------------------------------------------------------------------------

end

% formula based on matrix integration by parts
% Q' = expm(A*(h-s)), Q = -A^(-1)*expm(A*(h-s)) <- this might be missing a constant
% P = B*s, P' = B
% int(Q'*P,s,0,h)/h = (Q(h)*P(h) - Q(0)*P(0) - int(Q*P',s,0,h))/h
%  = (-A^(-1)*I*B*h - Q(0)*0 + A^(-1)*int(expm(A*(h-s))*B,s,0,h))/h
%  = -A^(-1)*B + A^(-1)/h*int(expm(A*(h-s))*B,s,0,h)
%  = A^(-1)*( int(expm(A*(h-s))*B,s,0,h)/h - B )