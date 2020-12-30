%--------------------------------------------------------------------------
% DTQP_reshape_X.m
% Reshape X into matrix form
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function X = DTQP_reshape_X(X,np,nt,ini)

% parameters
P = X(end-np+1:end);

% continuous variables into matrix form
X = reshape(X(1:end-np),nt,[]);

% combine with parameters, initial states, and final states in matrix form
X = [X,repelem(P',nt,1),repmat(X(1,ini{2}),nt,1),repmat(X(end,ini{2}),nt,1)];

end