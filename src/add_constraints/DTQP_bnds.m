%--------------------------------------------------------------------------
% DTQP_bnds.m
% Create sequences for simple upper and lower bound terms
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [I,V] = DTQP_bnds(bnd,p)
    % find time dependent matrix
    Bndt = DTQP_tmatrix(bnd.matrix,p);
    
    % permute
    Bndt = permute(Bndt,[1,3,2]);
    
    % variable locations for the variable type
    C = p.i{bnd.right};
    
    % initialize sequences
	I = []; V = [];
    
    % go through each variable of the current type
    for i = 1:length(C)
        
        switch bnd.right
            
            % control or states
            case {1,2} 
                
                % rows in lb/ub
                I = [I,DTQP_getQPIndex(C(i),bnd.right,1,p)];
                
                % nt values assigned
                V = [V,Bndt(i,:,:)];
                
            % parameters, initial states, or final states
            case {3,4,5}
                
                % row in lb/ub
                I = [I,DTQP_getQPIndex(C(i),bnd.right,0,p)];
                
                % single value assigned
                V = [V,Bndt(i,1,:)];
                
        end % end switch
        
    end % end for i
    
end % end function