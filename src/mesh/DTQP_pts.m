%--------------------------------------------------------------------------
% DTQP_pts.m
% Generate the time mesh (vector of discrete time values). Also, 
% potentially generate the quadrature weights and differentiation matrix
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [t,w,D] = DTQP_pts(in,dt)
% extract
t0 = in.t0; tf = in.tf;

% initialize
w = []; % empty
D = []; % empty

switch upper(dt.mesh)
    %----------------------------------------------------------------------
    case 'ED' % equidistant node points
        t = linspace(t0,tf,dt.nt)';
    %----------------------------------------------------------------------
    case 'LGL' % Legendre-Gauss-Lobatto nodes
        % scaled nodes and quadrature weights
        if strcmpi(dt.quadrature,'G') % Gaussian quadrature
            [tau,w] = lobpts(dt.nt); % using chebfun
        else % other
            tau = lobpts(dt.nt); % using chebfun
        end

        % differentiation matrix  
        if strcmpi(dt.defects,'PS')
            D = legslbdiff(dt.nt,tau);
        else
            D = []; % empty
        end

        % unscale mesh
        t = (tau + (tf+t0)/(tf-t0))*(tf-t0)*0.5;
    %----------------------------------------------------------------------
    case 'CGL' % Chebyshev-Gauss-Lobatto nodes
        % scaled nodes and quadrature weights
        if strcmpi(dt.quadrature,'CC') % Clenshaw-Curtis quadrature
            [tau,w] = chebpts(dt.nt); % using chebfun
        else % other
            tau = chebpts(dt.nt); % using chebfun
        end

        % differentiation matrix
        if strcmpi(dt.defects,'PS')
            D = diffmat(dt.nt,'chebkind2'); % using chebfun
        else
            D = []; % empty
        end

        % unscale mesh
        t = (tau + (tf+t0)/(tf-t0))*(tf-t0)*0.5;
    %----------------------------------------------------------------------
    case 'USER' % user-defined nodes
        if ~isfield(dt,'t')
            error('ERROR: opts.dt.t does not exist with USER option specified')
        else
            t = dt.t(:);
        end
        if strcmpi(dt.defects,'PS')
            error('PS option cannot handle USER mesh')
        end
    %----------------------------------------------------------------------
end

end