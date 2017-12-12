%--------------------------------------------------------------------------
% LGL_Dmatrix.m
% Determines approximate differentiation matrix for Legendre-based method
% with LGL nodes
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function D = DTQP_Dmatrix_LGL(tau)
    % number of nodes
    N = length(tau)-1;

    % See Page 110 of the book: J. Shen, T. Tang and L. Wang, Spectral Methods:
    % Algorithms, Analysis and Applications, Springer Series in Compuational
    % Mathematics, 41, Springer, 2011. 
    % Uses the function: lepoly()
    % Original function: D = legslbdiff(n,x) located at
    % http://www1.spms.ntu.edu.sg/~lilian/bookcodes/legen/legslbdiff.m
    n = N + 1;
    if n==0, D = []; return; end;   % null differentiation matrix
    xx = tau; y = lepoly(n-1,xx);
    D = (xx./y)*y'-(1./y)*(xx.*y)'; % compute L_{n-1}(x_j) (x_k-x_j)/L_{n-1}(x_k);     
                                    % 1/d_{kj} for k not= j (see (3.203)) 
    D = D + eye(n);                 % add the identity matrix so that 1./D can be operated                                     
    D = 1./D; 
    D = D - eye(n); 
    D(1,1) = -n*(n-1)/4; D(n,n) = -D(1,1); % update the diagonal entries  

end