%--------------------------------------------------------------------------
% DTQP_create_bnds.m
% Create simple lower and upper bound vectors
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [lb,ub] = DTQP_create_bnds(LB,UB,p)

    %----------------------------------------------------------------------
    % lower bounds
    %----------------------------------------------------------------------
    % initialize all lower bounds as -Inf (unconstrained)
    lb = -Inf*ones(p.nx,length(LB));
    
    % go through each substructure
    for i = 1:length(LB)
        
        % get the sequences for the lower bounds
        [I,V] = DTQP_bnds(LB(i),p);
        
        % combine
        lb(I,i) = V;
     
    end % end for i
    
    % assign maximum row value to lower bound vector
    lb = max(lb,[],2);

    %----------------------------------------------------------------------

    %----------------------------------------------------------------------
    % upper bounds
    %----------------------------------------------------------------------
    % initialize all upper bounds as Inf (unconstrained)
    ub = Inf*ones(p.nx,length(UB));

    % go through each substructure
    for i = 1:length(UB)
        
        % get the sequences for the upper bounds
        [I,V] = DTQP_bnds(UB(i),p);
        
        % combine
        ub(I,i) = V;
        
    end % end for i
    
    % assign minimum row value to lower bound vector
    ub = min(ub,[],2);
    %----------------------------------------------------------------------

end % end function