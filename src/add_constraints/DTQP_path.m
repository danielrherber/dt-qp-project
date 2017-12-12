%--------------------------------------------------------------------------
% DTQP_path.m
% Create sequences for path constraint terms
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [I,J,V,b] = DTQP_path(YZ,p)
    % initialize sequences
	I = []; J = []; V = [];
    
    % go through each substructure
    for j = 1:length(YZ.linear) % loop through the extended variables
        
        % find time dependent matrix
        YZt = DTQP_tmatrix(YZ.linear(j).matrix,p);
        
        % permute
        YZt = permute(YZt,[1,3,2]);
        
        % variable locations for the variable type
        C = p.i{YZ.linear(j).right};
        
        % row locations in 
        I = [I,repmat(1:p.nt,1,length(C))];
        
        % go through each variable of the current type
        for i = 1:length(C)
            
            % column locations
            J = [J,DTQP_getQPIndex(C(i),YZ.linear(j).right,1,p)]; 
            
            % nt values assigned
            V = [V,YZt(i,:,:)]; 
            
        end % end for i
        
    end % end for j
    
    % find time dependent vector
    bt = DTQP_tmatrix(YZ.b,p);

    % permute and assign
    b = permute(bt,[1,3,2]);

end % end function