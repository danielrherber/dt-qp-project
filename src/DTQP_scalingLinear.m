%--------------------------------------------------------------------------
% DTQP_scalingLinear.m
% Apply scaling to the problem based on a linear transformation
%--------------------------------------------------------------------------
% NOTE: constraint row scaling is currently commented
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [H,f,c,A,b,Aeq,beq,lb,ub,LLA,LRA,LLb,LRb,LLAeq,LRAeq,LLbeq,LRbeq,in,sm,sc] = ...
	DTQP_scalingLinear(H,f,c,A,b,Aeq,beq,lb,ub,LLA,LRA,LLb,LRb,LLAeq,LRAeq,LLbeq,LRbeq,in,s)

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
nr = length(sm);
Id = 1:nr;
sM = sparse(Id,Id,sm,nr,nr);

% scaling constant vector
sc = sparse([s1con;s2con;s3con]);

%--------------------------------------------------------------------------
% normal problem elements
%--------------------------------------------------------------------------
% c
if isempty(c)
	c = 0;
end

if ~isempty(find(f,1))
    c = c + f'*sc;
end

if ~isempty(find(H,1))
    c = c + sc'*H*sc/2;
end

% f
if ~isempty(find(f,1))
    f = sm.*f;
end

% H
if ~isempty(find(H,1))
    if ~isempty(find(f,1))
        f = f + sM*H*sc;
    else
        f = sM*H*sc;
    end

    H = sM*H*sM;
end

% A and b
if ~isempty(find(A,1))
    b = b - A*sc;
    A = A*sM;
end

% Aeq and beq
if ~isempty(find(Aeq,1))
    beq = beq - Aeq*sc;
    Aeq = Aeq*sM;
end

% lb
if ~isempty(find(lb,1))
    lb = (lb - sc)./sm;
end

% ub
if ~isempty(find(ub,1))
    ub = (ub - sc)./sm;
end

%--------------------------------------------------------------------------
% linkage constraints
%--------------------------------------------------------------------------
% LLA and LLb
if ~isempty(find(LLA,1))
    LLb = LLb - LLA*sc;
    LLA = LLA*sM;
end

% LRA and LRb
if ~isempty(find(LRA,1))
    LRb = LRb - LRA*sc;
    LRA = LRA*sM;
end

% LLAeq and LLbeq
if ~isempty(find(LLAeq,1))
    LLbeq = LLbeq - LLAeq*sc;
    LLAeq = LLAeq*sM;
end

% LRAeq and LRbeq
if ~isempty(find(LRAeq,1))
    LRbeq = LRbeq - LRAeq*sc;
    LRAeq = LRAeq*sM;
end

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