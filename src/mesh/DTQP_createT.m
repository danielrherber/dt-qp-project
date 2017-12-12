%--------------------------------------------------------------------------
% DTQP_createT.m
% Create time mesh (vector of discrete time values)
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function p = DTQP_createT(p,opts)

switch opts.NType
    %----------------------------------------------------------------------
    case 'ED' % equidistant node points
        p.t = linspace(p.t0,p.tf,p.nt)'; 
    %----------------------------------------------------------------------
    case 'LGL' % Lagrange-Gauss-Lobatto nodes
        tau = DTQP_nodes_LGL(p.nt-1);
        p.t = ( tau + (p.tf+p.t0)/(p.tf-p.t0) )*(p.tf-p.t0)/2;
    %----------------------------------------------------------------------
    case 'CGL' % Chebyshev-Gauss-Lobatto nodes
        tau = DTQP_nodes_CGL(p.nt-1);
        p.t = ( tau + (p.tf+p.t0)/(p.tf-p.t0) )*(p.tf-p.t0)/2;
    %----------------------------------------------------------------------
    case 'USER' % user-defined nodes
        if ~isfield(p,'t')
            error('ERROR: p.t does not exist with USER option specified')
        end
        if strcmp(opts.DTmethod,'PS')
            error('PS option cannot handle USER mesh')
        end
    %----------------------------------------------------------------------
end

end