%--------------------------------------------------------------------------
% DTQP_tmultiprod.m
% Evaluate the product of potentially time-varying matrices
%--------------------------------------------------------------------------
% - set 'prod' as the first entry to perform the matrix products
% - the matrices must be in the appropriate format for DTQP_tmatrix 
% - varargin can be used to provide a time mesh different than p.t
%--------------------------------------------------------------------------
% Example 1: two constant matrices ----------------------------------------
% A1 = rand(2,2); A2 = rand(2,2); p.t = (0:4)';
% matrices = {'prod',A1,A2};
% A = DTQP_tmultiprod(matrices,p);
% A(:,:,1) =
%     0.3754    0.2877
%     0.3754    0.2877
%     0.3754    0.2877
%     0.3754    0.2877
%     0.3754    0.2877
% A(:,:,2) =
%     0.3230    0.2350
%     0.3230    0.2350
%     0.3230    0.2350
%     0.3230    0.2350
%     0.3230    0.2350
% Example 2: three time-varying matrices ----------------------------------
% A1 = {@(t) sin(t),0;1,0}; A2 = {@(t) exp(-t) + 1,1;0,0}; A3 = eye(2);
% matrices = {'prod',A1,A2,A3}; t = (0:100)';
% A = DTQP_tmultiprod(matrices,p,t);
% A(1:4,:,:)
% ans(:,:,1) =
%          0    2.0000
%     1.1510    1.3679
%     1.0324    1.1353
%     0.1481    1.0498
% ans(:,:,2) =
%          0    1.0000
%     0.8415    1.0000
%     0.9093    1.0000
%     0.1411    1.0000
% The product equals:
% A1*A2*A3 = [ sin(t)*(exp(-t) + 1), sin(t)]
%            [          exp(-t) + 1,      1]
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function A = DTQP_tmultiprod(matrices,p,varargin)

if isempty(matrices)
    A = []; % empty result

elseif strcmp(matrices(1),'prod')
    % remove 'prod' specifier
    matrices(1) = [];

    % initial matrix
    A = DTQP_tmatrix(matrices{1},p,varargin{:});

    % compute the matrix products
    for k = 2:length(matrices)
        B = DTQP_tmatrix(matrices{k},p,varargin{:});
        A = multiprod(A,B,[2 3],[2 3]);
    end

else
    % only a single matrix (no matrix products needed)
    A = DTQP_tmatrix(matrices,p,varargin{:});
end

end