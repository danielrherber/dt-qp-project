%--------------------------------------------------------------------------
% DTQP_convolution.m
% Compute the convolution integral
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function Q = DTQP_convolution(A,B,p,opts)

% extract
t = p.t; h = p.h; nt = opts.dt.nt;

% number of states
ny = size(A,1);

% initialize convolution output matrix
Q = zeros(nt-1,ny,size(B,2));

%--------------------------------------------------------------------------
if isempty(B)
    return    
%--------------------------------------------------------------------------
elseif iscell(A)
        error('A should not be a function with the zero-hold method')
%--------------------------------------------------------------------------
% check if B is a cell (time-varying function)
elseif iscell(B)
    for k = 1:nt-1
        t0 = t(k); tf = t(k+1);
        if ~any(A(:)) % check if any elements of A are nonzero
            Q(k,:,:) = integral(@(tau) eye(ny)*squeeze(DTQP_tmatrix(B,p,tau))',...
                t0,tf,'ArrayValued',true,'RelTol',0,'AbsTol',1e-10);
        else % general case
            Q(k,:,:) = integral(@(tau) expm(A*(tf-tau))*squeeze(DTQP_tmatrix(B,p,tau))',...
                t0,tf,'ArrayValued',true,'RelTol',0,'AbsTol',1e-10);
        end
    end
%--------------------------------------------------------------------------
% check if A is singular
elseif rank(full(A)) < ny % abs(det(A)) < 100*eps
    if strcmpi(opts.dt.mesh,'ED') % equidistant points
        sysd  = c2d(ss(A,B,[],[]),h(1));
        q(1,:,:) = sysd.b;
        Q = repmat(q,[nt-1,1,1]);
    else % general meshes
        for k = 1:nt-1
            sysd  = c2d(ss(A,B,[],[]),h(k));
            Q(k,:,:) = sysd.b;
        end
    end
%--------------------------------------------------------------------------
% solve the matrix exponential integral then multiply
else
    % http://math.stackexchange.com/questions/658276/integral-of-matrix-exponential
    if strcmpi(opts.dt.mesh,'ED') % equidistant points
        v = A\( expm(A*(h(1))) - eye(ny) );
        q(1,:,:) = v*B;
        Q = repmat(q,[nt-1,1,1]);
    else % general meshes
        for k = 1:nt-1 % NOTE: could be replaced with multiprod
            v = A\( expm(A*(h(k))) - eye(ny) );
            Q(k,:,:) = v*B;
        end
    end
end

end