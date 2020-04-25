%--------------------------------------------------------------------------
% DTQP_scaling.m
% Apply scaling to the problem (linear and constraint row)
%--------------------------------------------------------------------------
% NOTE: constraint row scaling is currently commented
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [H,f,c,A,b,Aeq,beq,lb,ub,LAL,LAR,Lb,LAeqL,LAeqR,Lbeq,in,sm,sc] = ...
    DTQP_scaling(H,f,c,A,b,Aeq,beq,lb,ub,LAL,LAR,Lb,LAeqL,LAeqR,Lbeq,in,s)

%--------------------------------------------------------------------------
% START: simple scaling
%--------------------------------------------------------------------------
% extract
T = in.t; nt = in.nt; nu = in.nu; ny = in.ny; np = in.np;

% initialize as unity matrix scaling and zero constant shift
s1mat = ones(nu*nt,1); s2mat = ones(ny*nt,1); s3mat = ones(np,1);
s1con = zeros(nu*nt,1); s2con = zeros(ny*nt,1); s3con = zeros(np,1);

% obtain scaling vectors
for k = 1:length(s)

    % extract
    mat = s(k).matrix;
    sc = s(k).constant;

    switch s(k).right
        % controls
        case 1
            s1mat = createScalingVector(mat,nu,T,nt);
            s1con = createScalingVector(sc,nu,T,nt);
        % states
        case 2
            s2mat = createScalingVector(mat,ny,T,nt);
            s2con = createScalingVector(sc,ny,T,nt);
        % parameters
        case 3
            s3mat = createScalingVector(mat,np,[],1);
            s3con = createScalingVector(sc,np,[],1);
        % otherwise
        otherwise
            error(' ')
    end
end

% combine
sm = [s1mat;s2mat;s3mat];

% scaling diagonal matrix
nx = length(sm);
Is = 1:nx;
sM = sparse(Is,Is,sm,nx,nx);

% scaling constant vector
sc = sparse([s1con;s2con;s3con]);

%--------------------------------------------------------------------------
% c
if isempty(c)
	c = 0
end

if ~isempty(f)
    c = c + f'*sc;
end

if ~isempty(H)
    c = c + sc'*H*sc/2;
end

% f
if ~isempty(f)
    f = sm.*f;
end

% H
if ~isempty(H)
    if ~isempty(f)
        f = f + sM*H*sc;
    else
        f = sM*H*sc;
    end

    H = sM*H*sM;
end

% b
if ~isempty(b)
    b = b - A*sc;
end

% A
if ~isempty(A)
    A = A*sM;
end

% beq
if ~isempty(beq)
    beq = beq - Aeq*sc;
end

% Aeq
if ~isempty(Aeq)
    Aeq = Aeq*sM;
end

% lb
if ~isempty(lb)
    lb = (lb - sc)./sm;
end

% ub
if ~isempty(ub)
    ub = (ub - sc)./sm;
end

%--------------------------------------------------------------------------
% Lb
if ~isempty(Lb)
    Lb = Lb - LAL*sc; % not verified
end

% LAL
if ~isempty(LAL)
    LAL = LAL*sM; % not verified
end

% LAR
if ~isempty(LAR)
    LAR = LAR*sM; % not verified
end

% Lbeq
if ~isempty(Lbeq)
    Lbeq = Lbeq - LAeqL*sc; % not verified
end

% LAeqL
if ~isempty(LAeqL)
    LAeqL = LAeqL*sM; % not verified
end

% LAeqR
if ~isempty(LAeqR)
    LAeqR = LAeqR*sM; % not verified
end

%--------------------------------------------------------------------------
% END: simple scaling
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% START: constraint row scaling
%--------------------------------------------------------------------------
% r = 1./full(max(abs([A,b]),[],2));
% req = 1./full(max(abs([Aeq,beq]),[],2));
%
% % b
% if ~isempty(b)
%     b = b.*r;
% end
%
% % A
% if ~isempty(A)
%     A = bsxfun(@times,A,r);
% end
%
% % beq
% if ~isempty(beq)
%     beq = beq.*req;
% end
%
% % Aeq
% if ~isempty(Aeq)
%     Aeq = bsxfun(@times,Aeq,req);
% end

%--------------------------------------------------------------------------
% END: constraint row scaling
%--------------------------------------------------------------------------
end
% create the scaling vector for the particular input
function Y = createScalingVector(y,ny,T,nt)

if isa(y,'function_handle')
    % evaluate time-varying function
    Y = y(T);

    % reshape time-based matrix to column vector
    Y = Y(:);

elseif numel(y) == ny
    % expand scalar scaling
    Y = repelem(y(:),nt,1);

elseif all(size(y) == [nt ny])
    % reshape time-based matrix to column vector
    Y = y(:);

else
    error('wrong size')
end

end