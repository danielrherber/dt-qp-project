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
    % initialize storage arrays
    Isav = {}; Jsav = {}; Vsav = {};

    % go through each substructure
    for j = 1:length(YZ.linear) % loop through the extended variables

        % find time dependent matrix
        YZt = DTQP_tmatrix(YZ.linear(j).matrix,p);

        % variable locations for the variable type
        C = p.i{YZ.linear(j).right};

        % row locations
        Ir = (1:p.nt)';
        
        % go through each variable of the current type
        for i = 1:length(C)
            % row locations
            Is = Ir;

            % column locations
            Js = DTQP_getQPIndex(C(i),YZ.linear(j).right,1,p);

            % nt values assigned
            if length(size(YZt))==3
                Vs = YZt(:,:,i);
            else
                Vs = YZt(:,i,:);
            end

            % remove zeros
            ZeroIndex = find(~Vs);
            Is(ZeroIndex) = []; Js(ZeroIndex) = []; Vs(ZeroIndex) = [];

            % combine 
            Isav{end+1} = Is; Jsav{end+1} = Js; Vsav{end+1} = Vs;  

        end % end for i

    end % end for j

    % combine
    I = vertcat(Isav{:});
    J = vertcat(Jsav{:});
    V = vertcat(Vsav{:});
    
    % find time dependent vector
    b = DTQP_tmatrix(YZ.b,p);

end % end function