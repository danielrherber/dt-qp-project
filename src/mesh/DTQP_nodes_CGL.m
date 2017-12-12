%--------------------------------------------------------------------------
% DTQP_nodes_CGL.m
% Determine the Chebyshev-Gauss-Lobatto (CGL) nodes
%--------------------------------------------------------------------------
% tau = CGL_nodes(N)
%   N: number of nodes minus 1, should be an integer greater than 0
% tau: CGL nodes
%--------------------------------------------------------------------------
% Examples:
% tau = CGL_nodes(1)
% -1     1
% tau = CGL_nodes(2)
% -1     0     1
% tau = CGL_nodes(3)
% -1  -0.5   0.5   1
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function tau = DTQP_nodes_CGL(N)
%     % calculate node locations
%     tau = -cos(pi*(0:N)/n); % symbolically to maintain precision
%     tau = double(subs(tau,'n',N))';

    % calculate node locations
    tau = -cos(pi*(0:N)/N)';
end