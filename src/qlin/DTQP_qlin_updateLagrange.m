%--------------------------------------------------------------------------
% DTQP_qlin_updateLagrange.m
% Update Lagrange terms in the quasilinearization process
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function setup = DTQP_qlin_updateLagrange(setup,H,G,C,o,T,X,param)

% extract
nu = o.nu; ny = o.ny; np = o.np;

% convert to DTQP compatible functions
Hi = DTQP_qlin_update4tmatrix(H,T,X,param);
Gi = DTQP_qlin_update4tmatrix(G,T,X,param);
Ci = DTQP_qlin_update4tmatrix(C,T,X,param);

% initialize
keep = [];

% determine if nonzero entries are present (all variable types)
for k = 1:3
    I = DTQP_getProbIndex(k,nu,ny,np);
    if isnumeric(Hi)
        if any(Hi(I,:),'all') || any(Hi(:,I),'all')
            keep(end+1) = k;
        end
    else
        if any(~cellfun(@(x) isequal(x,0),Hi(I,:)),'all') || any(~cellfun(@(x) isequal(x,0),Hi(:,I)),'all')
            keep(end+1) = k;
        end
    end
    if isnumeric(Gi)
        if any(Gi(I),'all')
            keep(end+1:end+2) = [0,k];
        end
    else
        if any(~cellfun(@(x) isequal(x,0),Gi(I)),'all')
            keep(end+1:end+2) = [0,k];
        end
    end
end

% (TODO) check for constant term
% if any(Ci)
%     keep(end+1) = 0;
% end

% only need unique indices
keep = unique(keep);

% generate all permutations
[Iright,Ileft] = meshgrid(keep,keep);

% only consider when right is greater or equal to left
Ir = Iright(:) >= Ileft(:);
Iright = Iright(Ir);
Ileft = Ileft(Ir);

% initialize
L = struct('left',{},'right',{},'matrix',{});

% go through all entries
for k = 1:numel(Iright)

    % current matrix
    if (Ileft(k) ~= 0) && (Iright(k) ~= 0) % Hi
        M = Hi(DTQP_getProbIndex(Ileft(k),nu,ny,np),DTQP_getProbIndex(Iright(k),nu,ny,np));
    elseif (Ileft(k) == 0) && (Iright(k) == 0) % Ci
        M = Ci(1,1);
    else % Gi
        M = Gi(1,DTQP_getProbIndex(Iright(k),nu,ny,np));
    end

    % check if we want to add this term
    if isnumeric(M)
        if any(M,'all')
            L(end+1).matrix = M;
        else
            continue % don't add and continue
        end
    elseif any(~cellfun(@(x) isequal(x,0),M),'all')
        L(end+1).matrix = M;
    else
        continue % don't add and continue
    end

    % add left and right indices
    L(end).left = Ileft(k);
    L(end).right = Iright(k);
end

% add to setup structure
if isfield(setup,'L')
    setup.L = horzcat(setup.L,L); % this should not overwrite already included L
else
    setup.L = L;
end

end

% get optimization variables indices in original problem
function I = DTQP_getProbIndex(x,nu,ny,np)

switch x
    case 0 % singleton
        I = 1;
    case 1 % controls
        I = 1:nu;
    case 2 % states
        I = nu+1:nu+ny;
    case 3 % parameters
        I = nu+ny+1:nu+ny+np;
end

end