%--------------------------------------------------------------------------
% DTQP_reorder.m
% Reorder the optimization variables and constraint rows
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [H,f,c,A,b,Aeq,beq,lb,ub,p,X] = DTQP_reorder(H,f,c,A,b,Aeq,beq,lb,ub,p,X,flag)

if flag == 0
    %----------------------------------------------------------------------
    % START: reorder opt variables, [U,Y,p] -> [u_1,y_1,...,u_n,y_n,p]
    %----------------------------------------------------------------------
    % final index for continuous variables
    e = (p.nu+p.ns)*p.nt;
    
    % reshape to get sorting vector
    sV = reshape(reshape(1:e,[],p.nu+p.ns)',[],1);

    % add the parameters
    sV = [sV;e+1:e+p.np];
    
    % H
    if ~isempty(H)
        H = H(:,sV);
        H = H(sV,:);
    end
    
    % f
    if ~isempty(f)
        f = f(sV);
    end
    
    % A
    if ~isempty(A)
        A = A(:,sV);
    end

    % Aeq
    if ~isempty(Aeq)
        Aeq = Aeq(:,sV);
    end
    
    % lb
    if ~isempty(lb)
        lb = lb(sV);
    end

    % ub
    if ~isempty(ub)
        ub = ub(sV);
    end
    %----------------------------------------------------------------------
    % END: reorder opt variables, [U,Y,p] -> [u_1,y_1,...,u_n,y_n,p]
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % START: reorder linear constraint rows
    %----------------------------------------------------------------------
    if ~isempty(Aeq)    
        % sort the rows based on first nonzero entry
        [~,D] = sortrows(-abs(Aeq));
        
        % sort the matrix rows
        Aeq = Aeq(D,:);
        
        if ~isempty(beq)
            beq = beq(D);
        end
    end

    if ~isempty(A)    
        % sort the rows based on first nonzero entry
        [~,D] = sortrows(-abs(A));
        
        % sort the matrix rows
        A = A(D,:);
        
        if ~isempty(b)
            b = b(D);
        end
    end
    %----------------------------------------------------------------------
    % END: reorder linear constraint rows
    %----------------------------------------------------------------------

else
    % original ordering
    % reorder, [u_1, y_1, ... , u_n, y_n, p] -> [U,Y,p]
    sV = [];
    for i = 1:(p.nu+p.ns)
        sV = [sV, i:(p.nu+p.ns):p.nx];
    end

    % unsort
    X = X(sV);

end

end