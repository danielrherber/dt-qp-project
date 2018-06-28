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
function [J,V,b] = DTQP_boundary(yz,in)
    % extract some of the variables
    nt = in.nt; nu = in.nu; ny = in.ny; ini = in.i;

    % initialize storage arrays
    Jsav = cell(0,1); Vsav = cell(0,1);

    % go through each substructure
    for j = 1:length(yz.linear) % loop through the extended variables

        % check if the supplied matrix is a cell (column vector)
        if ~ismatrix(yz.linear(j).matrix)
            yzt = cell2mat(yz.linear(j).matrix);
        else
            yzt = yz.linear(j).matrix;
        end

        % variable locations for the variable type
        C = ini{yz.linear(j).right};

        % go through each variable of the current type
        for i = 1:length(C)

            % column location
            Js = DTQP_getQPIndex(C(i),yz.linear(j).right,0,nt,nu,ny);

            % single value assigned
            Vs = reshape(yzt(i,1),1,[]);

            % remove zeros
            ZeroIndex = find(~Vs);
            Js(ZeroIndex,:) = []; Vs(ZeroIndex,:) = [];

            % combine 
            Jsav{end+1} = Js; Vsav{end+1} = Vs;  

        end % end for i

    end % end for j

    % combine
    J = vertcat(Jsav{:});
    V = vertcat(Vsav{:});
    
    % assign constant value
    b = yz.b;

end % end function