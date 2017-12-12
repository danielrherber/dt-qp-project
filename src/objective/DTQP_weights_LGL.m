%--------------------------------------------------------------------------
% DTQP_weights_LGL.m
% Determines Gaussian quadrature weights using Lagrange-Gauss-Lobatto (LGL)
% nodes
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function w = DTQP_weights_LGL(tau)
    % number of nodes
    N = length(tau)-1;
    
    % See Page 99 of the book: J. Shen, T. Tang and L. Wang, Spectral Methods:
    % Algorithms, Analysis and Applications, Springer Series in Compuational
    % Mathematics, 41, Springer, 2011. 
    % Uses the function: lepoly() 
    % Original function: [varargout] = legslb(n) located at
    % http://www1.spms.ntu.edu.sg/~lilian/bookcodes/legen/legslb.m
    [~,y] = lepoly(N,tau(2:end-1));
    % Use the weight expression (3.188) to compute the weights
    w = [2/(N*(N+1));2./(N*(N+1)*y.^2);2/(N*(N+1))]; 

end