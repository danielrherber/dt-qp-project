%--------------------------------------------------------------------------
% DTQP_boundary.m
% Create sequences for boundary constraint terms
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [J,V,b] = DTQP_boundary(yz,in)

% extract some of the variables
nt = in.nt; ini = in.i; I_stored = in.I_stored;

% initialize storage arrays
Jsav = zeros(0,1); Vsav = zeros(0,1);

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

        % single value assigned
        Vs = yzt(i,1);

        % check if any entries are nonzero
        if any(Vs)

            % column location
            Js = DTQP_getQPIndex(C(i),yz.linear(j).right,0,nt,I_stored);

            % combine
            Jsav(end+1,1) = Js; Vsav(end+1,1) = Vs;

        end

    end % end for i

end % end for j

% combine
J = Jsav;
V = Vsav;

% assign constant value
b = yz.b;

end % end function