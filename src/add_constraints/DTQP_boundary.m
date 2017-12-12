%--------------------------------------------------------------------------
% DTQP_boundary.m
% Create sequences for boundary constraint terms
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [J,V,b] = DTQP_boundary(yz,p)
    % initialize sequences
	J = []; V = [];

    % go through each substructure
    for j = 1:length(yz.linear) % loop through the extended variables

        % check if the supplied matrix is a cell
        if ~ismatrix(yz.linear(j).matrix)
            yzt = cell2mat(yz.linear(j).matrix);
        else
            yzt = yz.linear(j).matrix;
        end

        % variable locations for the variable type
        C = p.i{yz.linear(j).right};

        % go through each variable of the current type
        for i = 1:length(C)

                % column location
                J = [J,DTQP_getQPIndex(C(i),yz.linear(j).right,0,p)];

                % single value assigned
                V = [V,reshape(yzt(i,1),1,[])];
                
        end % end for i
        
    end % end for j

    % assign constant value
    b = yz.b;

end % end function