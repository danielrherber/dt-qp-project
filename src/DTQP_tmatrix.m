%--------------------------------------------------------------------------
% DTQP_tmatrix.m
% Evaluate a potentially time-varying matrix on a specific time mesh 
%--------------------------------------------------------------------------
% - the matrix must be vectorized if time-varying
% - varargin can be used to provide a time mesh different than p.t
%--------------------------------------------------------------------------
% Example 1: basic time-invariant -----------------------------------------
% p.t = 0:4; At = DTQP_tmatrix(ones(2),p);
% At(:,:,1) =
%      1     1
%      1     1
%      1     1
%      1     1
%      1     1
% At(:,:,2) =
%      1     1
%      1     1
%      1     1
%      1     1
%      1     1
% Example 2: basic time-varying -------------------------------------------
% p.t = 0:2; A{1} = @(t) sin(t); A{2} = -1; At = DTQP_tmatrix(A,p);
% At(:,:,1) =
%          0
%     0.8415
%     0.9093
% At(:,:,2) =
%     -1
%     -1
%     -1
% Example 3: passing additional parameters --------------------------------
% p.t = 0:2; p.a = 3; A{1} = @(t,p) p.a*sin(t); At = DTQP_tmatrix(A,p);
% At =
%          0
%     2.5244
%     2.7279
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function At = DTQP_tmatrix(A,p,varargin)

    % check if another time mesh is inputted
    if isempty(varargin)
        t = p.t;
    else
        t = varargin{1};
    end

    if isempty(A) % A is empty
        At = []; % return empty matrix if A is empty
    elseif isnumeric(A) % constant matrix
        % add dummy dimension
        A1(1,:,:) = full(A); % needs to be full
        
        % replicate the matrix nt times
        At = repmat(A1,[length(t),1,1]);

    else

        % convert to a cell if a normal matrix or single time-varying function
        if isa(A,'function_handle')
            A = {A};
        elseif ~isa(A,'cell')
            A = num2cell(A);
        elseif isequal(size(A),[1 1]) && ~isa(A{1},'function_handle')
            A = num2cell(A{1});
        end

        % get matrix size
        r = size(A,1);
        c = size(A,2);

        % initialize output matrix
        At = zeros(length(t),r,c);

        % go through each row and column in A
        for i = 1:r
            for j = 1:c
                if isempty(A{i,j}) % A is empty
                    % do nothing
                elseif isa(A{i,j},'function_handle') % A is time-varying
                    if nargin(A{i,j}) == 2
                        At(:,i,j) = A{i,j}(t,p);
                    elseif nargin(A{i,j}) == 1
                        At(:,i,j) = A{i,j}(t);
                    else
                        error('not a properly defined function')
                    end
                elseif A{i,j} == 0 % A is zero
                    % do nothing
                else % A is time-invariant
                    At(:,i,j) = A{i,j};
                end
            end   
        end
    end
end