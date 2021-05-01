%--------------------------------------------------------------------------
% DTQP_reorder.m
% Reorder the optimization variables and constraint rows
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = DTQP_reorder(in,varargin)

if length(varargin) == 9
    % extract
    H = varargin{1}; f = varargin{2}; c = varargin{3};
    A = varargin{4}; b = varargin{5};
    Aeq = varargin{6}; beq = varargin{7};
    lb = varargin{8}; ub = varargin{9};

    %----------------------------------------------------------------------
    % START: reorder opt variables, [U,Y,p] -> [u_1,y_1,...,u_n,y_n,p]
    %----------------------------------------------------------------------
    % final index for continuous variables
    e = (in.nu+in.ny)*in.nt;

    % reshape to get sorting vector
    sV = reshape(reshape(1:e,[],in.nu+in.ny)',[],1);

    % add the parameters
    sV = [sV;e+1:e+in.np];

    % H
    if ~isempty(find(H,1))
        H = H(:,sV);
        H = H(sV,:);
    end

    % f
    if ~isempty(find(f,1))
        f = f(sV);
    end

    % A
    if ~isempty(find(A,1))
        A = A(:,sV);
    end

    % Aeq
    if ~isempty(find(Aeq,1))
        Aeq = Aeq(:,sV);
    end

    % lb
    if ~isempty(find(lb,1))
        lb = lb(sV);
    end

    % ub
    if ~isempty(find(ub,1))
        ub = ub(sV);
    end
    %----------------------------------------------------------------------
    % END: reorder opt variables, [U,Y,p] -> [u_1,y_1,...,u_n,y_n,p]
    %----------------------------------------------------------------------

    %----------------------------------------------------------------------
    % START: reorder linear constraint rows
    %----------------------------------------------------------------------
    if ~isempty(find(Aeq,1))

        % % sort the rows based on first nonzero entry
        % [~,D] = sortrows(-abs(Aeq));

        % number of constraints
        e = size(Aeq,1);

        % reshape to get sorting vector
        D = reshape(reshape(1:e,[],in.ny)',[],1);

        % sort the matrix rows
        Aeq = Aeq(D,:);

        if ~isempty(find(beq,1))
            beq = beq(D);
        end

    end

    if ~isempty(find(A,1))

        % % sort the rows based on first nonzero entry
        % [~,D] = sortrows(-abs(A));

        % number of constraints
        e = size(A,1);

        % reshape to get sorting vector
        D = reshape(reshape(1:e,[],in.ny)',[],1);

        % sort the matrix rows
        A = A(D,:);

        if ~isempty(find(b,1))
            b = b(D);
        end

    end
    %----------------------------------------------------------------------
    % END: reorder linear constraint rows
    %----------------------------------------------------------------------

    % assign outputs
    varargout{1} = H; varargout{2} = f; varargout{3} = c;
    varargout{4} = A; varargout{5} = b;
    varargout{6} = Aeq; varargout{7} = beq;
    varargout{8} = lb; varargout{9} = ub;

else

    % original ordering
    % reorder, [u_1, y_1, ... , u_n, y_n, p] -> [U,Y,p]
    sV = [];
    for i = 1:(in.nu+in.ny)
        sV = [sV, i:(in.nu+in.ny):in.nx];
    end
    % sV = reshape(1:in.nx,in.nu+in.ny,[])'; % alternative implementation

    % unsort
    varargout{1} = varargin{1}(sV);

end

end