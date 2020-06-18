%--------------------------------------------------------------------------
% DTQP_scalingRows.m
% Apply scaling to the constraint rows
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [A,b,Aeq,beq] = DTQP_scalingRows(A,b,Aeq,beq)

% determine row scaling factors
r = 1./max(abs([A,b]),[],2);
req = 1./max(abs([Aeq,beq]),[],2);

% construct sparse diagonal matrices
nr = size(A,1); % number of rows
Id = 1:nr; % indices of diagonal
rM = sparse(Id,Id,r,nr,nr);

nr = size(Aeq,1); % number of rows
Id = 1:nr; % indices of diagonal
reqM = sparse(Id,Id,req,nr,nr);

% b
if ~isempty(find(b,1))
    b = r.*b;
end

% A
if ~isempty(find(A,1))
    A = rM*A;
end

% beq
if ~isempty(find(beq,1))
    beq = req.*beq;
end

% Aeq
if ~isempty(find(Aeq,1))
    Aeq = reqM*Aeq;
end

end