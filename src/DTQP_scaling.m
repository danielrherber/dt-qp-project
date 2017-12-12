%--------------------------------------------------------------------------
% DTQP_scaling.m
% Apply scaling to the problem (linear and constraint row)
%--------------------------------------------------------------------------
% NOTE: constraint row scaling is currently commented
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [H,f,A,b,Aeq,beq,lb,ub,p,s] = DTQP_scaling(H,f,A,b,Aeq,beq,lb,ub,p,s)
    
    %----------------------------------------------------------------------
    % START: simple scaling
    %----------------------------------------------------------------------
    % initialize as unity linear scaling
    s1 = ones(p.nu*p.nt,1); s2 = ones(p.ns*p.nt,1); s3 = ones(p.np,1);

    % scale controls
    for k = 1:length(s)
        % extract matrix
        m = s(k).matrix(:);
        
        switch s(k).right
            % controls
            case 1
                if p.nu == length(m)
                    s1 = kron(m,ones(p.nt,1));
                else
                    error('wrong size')
                end
            % states
            case 2
                if p.ns == length(m)
                    s2 = kron(m,ones(p.nt,1));
                else
                    error('wrong size')
                end
            % parameters
            case 3
                if p.np == length(m)
                    s3 = m;
                else
                    error('wrong size')
                end
            % otherwise
            otherwise
                error(' ')
        end
    end
    
    % combine
    s = [s1;s2;s3];
    
    % diagonal scaling matrix
    S = sparse(diag(s));    

    % scale the matrices
    % H
    if ~isempty(H)
        H = S*H*S;
    end
    
    % f
    if ~isempty(f)
        f = s.*f;
    end
    
    % A
    if ~isempty(A)
        A = A*S;
    end

    % Aeq
    if ~isempty(Aeq)
        Aeq = Aeq*S;
    end
    
    % lb
    if ~isempty(lb)
        lb = lb./s;
    end 
    
    % ub
    if ~isempty(ub)
        ub = ub./s;
    end 
    
    % Note: b and beq are not scaled  
    
    %----------------------------------------------------------------------
    % END: simple scaling
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % START: constraint row scaling
    %----------------------------------------------------------------------
%     r = max(abs(A),[],2);
%     req = max(abs(Aeq),[],2);
%     
%     b
%     if ~isempty(b)
%         b = b./r;
%     end
%     
%     A
%     if ~isempty(A)
%         A = bsxfun(@rdivide,A,r);
%     end
%     
%     beq
%     if ~isempty(beq)
%         beq = beq./req;
%     end
%     
%     Aeq
%     if ~isempty(Aeq)
%         Aeq = bsxfun(@rdivide,Aeq,req);
%     end

    %----------------------------------------------------------------------
    % END: constraint row scaling
    %----------------------------------------------------------------------
end